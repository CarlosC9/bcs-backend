

# ALGORITHMS
import datetime

from sqlalchemy import Integer, Column, String, Text, ForeignKey, JSON, DateTime
from sqlalchemy.orm import relationship, backref

from biobarcoding.db_models import ORMBase, GUID

prefix = "jobs_"


class AlgorithmType(ORMBase):
    """ """
    __tablename__ = f"{prefix}algorithm_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class Algorithm(ORMBase):
    """ Description of a specific algorithm"""
    __tablename__ = f"{prefix}algorithms"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    description = Column(Text)


class Workflow(ORMBase):
    """ Bioinformatic repeatable workflows executed in a compute resource to process bioinformatic information """
    __versioned__ = {}
    __tablename__ = f"{prefix}wfs"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    description = Column(Text)
    schema_inputs = Column(JSON)  # How to specify inputs for the workflow
    schema_outputs = Column(JSON)  # How to read workflow outputs
    execution = Column(JSON)  # How to call the workflow in a resource independent ("generic") way


class WorkflowAlgorithm(ORMBase):
    """ One of the steps in a Workflow """
    __tablename__ = f"{prefix}wfs_algorithms"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    workflow_id = Column(Integer, ForeignKey(Workflow.id), nullable=False, primary_key=False)
    workflow = relationship(Workflow, backref=backref("algorithms", cascade="all, delete-orphan"))


class WorkflowStep(ORMBase):
    """ One of the steps in a Workflow """
    __tablename__ = f"{prefix}wfs_steps"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    workflow_id = Column(Integer, ForeignKey(Workflow.id), nullable=False, primary_key=False)
    workflow = relationship(Workflow, backref=backref("steps", cascade="all, delete-orphan"))


class ProcessorArchitecture(ORMBase):
    """ Processor and architecture """
    __tablename__ = f"{prefix}processor_architectures"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class OperatingSystem(ORMBase):
    """ Operating system code """
    __tablename__ = f"{prefix}operating_systems"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class ComputingType(ORMBase):
    """ Computing types: sequencial, MPI, OpenMP, GPU, ... """
    __tablename__ = f"{prefix}computing_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class JobManagementType(ORMBase):
    __tablename__ = f"{prefix}job_mgmt_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class ComputeResource(ORMBase):
    """ Specific compute resource """

    __tablename__ = f"{prefix}compute_resources"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    proc_arch_id = Column(Integer, ForeignKey(ProcessorArchitecture.id), nullable=False, primary_key=False)
    proc_arch = relationship(ProcessorArchitecture)
    os_id = Column(Integer, ForeignKey(OperatingSystem.id), nullable=False, primary_key=False)
    operating_system = relationship(OperatingSystem)
    jm_type_id = Column(Integer, ForeignKey(JobManagementType.id), nullable=False, primary_key=False)
    jm_type = relationship(JobManagementType)
    jm_credentials = Column(JSON)
    agreement = Column(JSON)
    consumption_counters = Column(JSON)


class WorkflowInComputeResource(ORMBase):
    """
    Specific adaptation of a Worflow to a ComputeResource
    """

    __tablename__ = f"{prefix}wfs_compute_resources"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    workflow_id = Column(Integer, ForeignKey(Workflow.id), nullable=False, primary_key=False)
    workflow = relationship(Workflow, backref=backref("resource_adaptations", cascade="all, delete-orphan"))
    resource_id = Column(Integer, ForeignKey(ComputeResource.id), nullable=False, primary_key=False)
    resource = relationship(ComputeResource, backref=backref("workflow_adaptations", cascade="all, delete-orphan"))


class JobStatus(ORMBase):
    """ created, preparing_workspace, transferring_data_to_resource, submitted, completed_successfully, transferring_data_from_resource, cleaning_up_workspace,
    cancelled, paused, completed_error """
    __tablename__ = f"{prefix}job_statuses"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class Job(ORMBase):
    """ Core entity in the ALgorithmic domain.
    Workflows are launched into ComputeResources with some parameters and under an Identity
    Then, jobs can be checked, cancelled and results obtained
    Also, listing of filtered/ordered jobs may be retrieved
    """
    __tablename__ = f"{prefix}jobs"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)

    workflow_id = Column(Integer, ForeignKey(Workflow.id), nullable=False, primary_key=False)
    workflow = relationship(Workflow)
    resource_id = Column(Integer, ForeignKey(ComputeResource.id), nullable=False, primary_key=False)
    resource = relationship(ComputeResource)
    log = Column(Text)
    inputs = Column(JSON)
    outputs = Column(JSON)


class JobStatusLog(ORMBase):
    __tablename__ = f"{prefix}jobs_status_log"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)

    job_id = Column(Integer, ForeignKey(Job.id), nullable=False, primary_key=False)
    job = relationship(Job)
    logged_time = Column(DateTime, default=datetime.datetime.utcnow())
    logged_status_id = Column(Integer, ForeignKey(JobStatus.id), nullable=False, primary_key=False)
    logged_status = relationship(JobStatus, backref=backref("statuses", cascade="all, delete-orphan"))
