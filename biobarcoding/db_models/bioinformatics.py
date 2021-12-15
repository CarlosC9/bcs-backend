# BIOINFORMATIC data
import uuid

from sqlalchemy import Column, Integer, ForeignKey, UnicodeText, String, BigInteger
from sqlalchemy.orm import relationship

from . import ORMBase, GUID
from .core import data_object_type_id, FunctionalObject

prefix = "bo_"

data_object_type_id.update({
    "sequence": 11,
    "multiple-sequence-alignment": 12,
    "phylogenetic-tree": 13,
    "sequence-similarity": 14,  # BLAST
})


class BarCodingRegions(ORMBase):
    __tablename__ = f"{prefix}bar_coding_regions"
    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    description = Column(UnicodeText)


class Specimen(ORMBase):
    __tablename__ = f"{prefix}specimens"
    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    description = Column(UnicodeText)


class Sequence(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}sequences"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['sequence'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)
    specimen_id = Column(Integer, ForeignKey(Specimen.id), nullable=True, primary_key=False)
    specimen = relationship(Specimen)


class MultipleSequenceAlignment(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}msas"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['multiple-sequence-alignment'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)


class PhylogeneticTree(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}phylo_trees"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['phylogenetic-tree'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)


class SequenceSimilarity(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}seq_sims"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['sequence-similarity'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)
