# BIOINFORMATIC data
from sqlalchemy import Column, Integer, ForeignKey, UnicodeText, Unicode, String

from biobarcoding.db_models import ORMBase, GUID, ObjectType

prefix = "bo_"


class BioinformaticObject(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}bos"
    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    bo_type_id = Column(Integer, ForeignKey(ObjectType.id))
    name = Column(String(80))
    content = Column(UnicodeText)

    __mapper_args__ = {
        'polymorphic_identity': 'bioinformatic_obj',
        'polymorphic_on': bo_type_id
    }


class BarCodingRegions(ORMBase):
    __tablename__ = f"{prefix}bar_coding_regions"
    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    description = Column(UnicodeText)


class Sequence(BioinformaticObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}sequences"
    __mapper_args__ = {
        'polymorphic_identity': 'sequence',
    }

    id = Column(Integer, ForeignKey(BioinformaticObject.id), primary_key=True)


class MultipleSequenceAlignment(BioinformaticObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}msas"
    __mapper_args__ = {
        'polymorphic_identity': 'msa',
    }

    id = Column(Integer, ForeignKey(BioinformaticObject.id), primary_key=True)


class PhylogeneticTree(BioinformaticObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}phylo_trees"
    __mapper_args__ = {
        'polymorphic_identity': 'phylo-tree',
    }

    id = Column(Integer, ForeignKey(BioinformaticObject.id), primary_key=True)
