import datetime
import re
import uuid

from sqlalchemy import Column, Integer, ForeignKey, String, BigInteger, Boolean, DateTime, UniqueConstraint, event, \
    func, Index, or_, and_
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship, backref
from sqlalchemy_utils import TSVectorType

from . import ORMBase, GUID, ObjectType
from .jobs import Job
from .. import app_acronym
from ..common import generate_json

prefix = f"{app_acronym}_"

data_object_type_id = {
    "dataset": 0,
    "geolayer": 1,
    "dataframe": 2,
    "grid": 3,
    "process": 100,
    "geoprocess": 101,
    "geoprocess_instance": 102,
    "case_study": 103
}


class FunctionalObject(ORMBase):
    __tablename__ = f"{prefix}fos"
    id = Column(BigInteger, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)

    obj_type_id = Column("do_type_id", Integer, ForeignKey(ObjectType.id))
    native_id = Column(BigInteger, nullable=True, primary_key=False)
    native_table = Column(String(80))

    name = Column(String(300))
    attributes = Column(JSONB)  # Tags, and other categories, used to classify the object

    ts_vector = Column(TSVectorType())
    ts_vector_update_time = Column(DateTime, default=datetime.datetime.utcnow())
    entity_update_time = Column(DateTime, default=datetime.datetime.utcnow())

    __table_args__ = (
        UniqueConstraint(native_table, native_id, name=__tablename__ + '_c1'),
        Index('ix_fobj__ts_vector__', ts_vector, postgresql_using='gin'),
    )

    __mapper_args__ = {
        'polymorphic_identity': 'functional_obj',
        'polymorphic_on': obj_type_id
    }


def set_functional_object_tsvector(entity):
    """
    Prepare and set the tsvector for the entity

    :param entity:
    :return:
    """
    # Name
    if entity.name:
        lst_parts = re.split(r"\s|(?<!\d)[,.](?!\d)", entity.name)
    else:
        lst_parts = []
    # ID
    lst_parts.append(str(entity.id))
    # Attributes
    if entity.attributes is not None:
        lst_parts.append(generate_json(entity.attributes))
    if isinstance(entity, Dataset):
        # Creation time
        lst_parts.append(str(entity.creation_time))
        if entity.structure is not None:
            lst_parts.append(generate_json(entity.structure))
        if entity.provenance is not None:
            lst_parts.append(generate_json(entity.provenance))
        from ..db_models.geographics import GeographicLayer
        if isinstance(entity, GeographicLayer):
            if entity.properties is not None:
                lst_parts.append(generate_json(entity.properties))
            if entity.layer_type is not None:
                lst_parts.append(str(entity.layer_type))

    entity.ts_vector = func.to_tsvector(' '.join([i for i in lst_parts if i is not None]))  # !!
    entity.ts_vector_update_time = datetime.datetime.utcnow()


def after_create_or_update(mapper, connection, target):
    if target.entity_update_time is None or \
            (target.entity_update_time is not None and
             (target.ts_vector_update_time is None or
              target.ts_vector_update_time > target.entity_update_time)):
        set_functional_object_tsvector(target)


event.listen(FunctionalObject, 'before_insert', after_create_or_update, propagate=True)
event.listen(FunctionalObject, 'before_update', after_create_or_update, propagate=True)


# Refresh all!!
def update_functional_object_tsvector(session):
    qry = session.query(FunctionalObject).\
        filter(or_(FunctionalObject.entity_update_time == None,
                   and_(FunctionalObject.entity_update_time != None,
                        or_(FunctionalObject.ts_vector_update_time == None,
                            FunctionalObject.ts_vector_update_time < FunctionalObject.entity_update_time)))).all()
    for entity in qry:
        entity.entity_update_time = datetime.datetime.utcnow()
        set_functional_object_tsvector(entity)
    session.commit()


class CaseStudy(FunctionalObject):
    """
    A case study is a set of Functional Objects, with the purpose of reducing the size of the list of Functional Objects
    when browsing. It may also have permissions attached to it

    """
    __tablename__ = f"{prefix}case_studies"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['case_study'],
    }

    id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)


class CaseStudy2FunctionalObject(ORMBase):
    """ Many to Many table, to relate Case Studies and FunctionalObjects """

    __tablename__ = f"{prefix}case_studies_functional_objects"

    case_study_id = Column(BigInteger, ForeignKey(CaseStudy.id), primary_key=True)
    functional_object_id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)
    case_study = relationship(CaseStudy, foreign_keys=[case_study_id], backref=backref("fos", cascade="all, delete-orphan"))
    functional_object = relationship(FunctionalObject, foreign_keys=[functional_object_id], backref=backref("case_studies", cascade="all, delete-orphan"))


class Dataset(FunctionalObject):
    """
    Root of the dataset hierarchy of objects, whose only member currently will be "GeographicLayer"
    Other types like Grid or DataFrame may be incorporated in future improvements
    """
    __tablename__ = f"{prefix}datasets"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['dataset'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)
    identity_id = Column(Integer)
    creation_time = Column(DateTime, default=datetime.datetime.utcnow())
    structure = Column(JSONB)  # Used to match datasets and ports
    provenance = Column(JSONB)  # Origin of the dataset
    # attributes: source/intermediate/output (source+output, intermediate+output are possible), provenance, ...


class PortType(ORMBase):
    """ A Port Type is used to match Layers with CProcess ports
        If GeographicLayers and CProcessPort "PortType" are the same, they match
    """
    __tablename__ = f"{prefix}port_types"
    id = Column(BigInteger, primary_key=True, autoincrement=True)
    name = Column(String(128), nullable=True)
    # Fields characterizing the port type
    attributes = Column(JSONB)


class DatasetTypes(ORMBase):
    """ Many to Many table, to relate Datasets and PortTypes """

    __tablename__ = f"{prefix}dataset_port_types"

    dataset_id = Column(BigInteger, ForeignKey(Dataset.id), primary_key=True)
    port_type_id = Column(BigInteger, ForeignKey(PortType.id), primary_key=True)
    dataset = relationship(Dataset, backref=backref("types", cascade="all, delete-orphan"))
    port_type = relationship(PortType)


class CProcess(FunctionalObject):
    """
    GeoProcess Definition
    """
    __tablename__ = f"{prefix}processes"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['geoprocess'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)
    # Information on how to execute
    process_type = Column(String(80))  # PostGIS, Python, R, ...


class CProcessPort(ORMBase):
    """ A port in a GeoProcess Definition """
    __tablename__ = f"{prefix}processes_ports"
    id = Column(BigInteger, primary_key=True, autoincrement=True)

    name = Column(String(128))
    process_id = Column(BigInteger, ForeignKey(CProcess.id))
    process = relationship(CProcess, backref=backref("ports", cascade="all, delete-orphan"))
    # To unify, ports of the same type point to a "PortType" register
    port_type_id = Column(BigInteger, ForeignKey(PortType.id), nullable=True)
    port_type = relationship(PortType)
    input = Column(Boolean, nullable=False)  # True: Input; False: Output
    # Other fields characterizing the port
    attributes = Column(JSONB)


class CProcessInstance(FunctionalObject):
    """ Point to Process.
    Can be "candidate for execution", "in progress" or finished (successfully or not)
    Different from Job which is in charge of carrying out the process itself
    """
    __versioned__ = {}
    __tablename__ = f"{prefix}process_instances"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['geoprocess_instance'],
    }

    id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)

    instantiated_process_id = Column(BigInteger, ForeignKey(CProcess.id))
    instantiated_process = relationship(CProcess, foreign_keys=[instantiated_process_id])
    #
    # STATUSES: "candidate",
    #           "scheduled" (has "job_id", the specific status is delegated to Job class)
    #           "cancelled", (Job execution cancelled, do not set directly)
    #           "success" (completed successfully),
    #           "error" (completed unsuccessfully)
    #
    status = Column(String(10))
    # If it has a value, the process has been launched for execution
    job_id = Column(BigInteger, ForeignKey(Job.id))
    creation_time = Column(DateTime, default=datetime.datetime.utcnow())
    hash_of_canonical = Column(String(64), unique=True)
    params = Column(JSONB)


class PortInProcessInstance(ORMBase):
    """ Point to ProcessExecution, to Dataset and to Port """
    __tablename__ = f"{prefix}process_instances_ports"
    process_instance_id = Column(BigInteger, ForeignKey(CProcessInstance.id), primary_key=True)
    process_instance = relationship(CProcessInstance, backref=backref("ports", cascade="all, delete-orphan"))
    port_id = Column(BigInteger, ForeignKey(CProcessPort.id), primary_key=True)
    port = relationship(CProcessPort)
    dataset_id = Column(BigInteger, ForeignKey(Dataset.id), nullable=True)
    dataset = relationship(Dataset)
