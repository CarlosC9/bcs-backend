# BIOINFORMATIC data
import uuid

from sqlalchemy import Column, Integer, ForeignKey, UnicodeText, String, BigInteger, UniqueConstraint
from sqlalchemy.orm import relationship, backref

from . import ORMBase, GUID
from .core import data_object_type_id, FunctionalObject
from .metadata import Taxon

prefix = "bo_"

analysis_object_type_id = {
    "multiple-sequence-alignment": 12,
    "phylogenetic-tree": 13,
    "sequence-similarity": 14,  # BLAST
    "discriminant-matrix": 15,
    "supermatrix": 16,
}
data_object_type_id.update({
    "specimen": 10,
    "sequence": 11,
    "analysis": list(analysis_object_type_id.values()),
    **analysis_object_type_id
})


class BarCodingRegions(ORMBase):
    __tablename__ = f"{prefix}bar_coding_regions"
    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    description = Column(UnicodeText)


class Gene(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}genes"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    description = Column(UnicodeText)


class GeneSynonym(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}gene_synonyms"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    name = Column(String(80))
    gene_id = Column(BigInteger, ForeignKey(Gene.id), nullable=False)
    gene = relationship(Gene, backref=backref("synonyms", cascade="all, delete-orphan"))


class Locus(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}loci"

    # TODO:
    #  Cytogenetic location <chromosome/arm/region/band> p.e. 1q1.2-q2.3
    #  Molecular location, nucleotides positions
    id = Column(BigInteger, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    left = Column(Integer)
    right = Column(Integer)
    gene_id = Column(BigInteger, ForeignKey(Gene.id), nullable=True)
    gene = relationship(Gene, backref=backref("locus", cascade="all, delete-orphan"))
    taxon_id = Column(BigInteger, ForeignKey(Taxon.id), nullable=False)
    taxon = relationship(Taxon, backref=backref("locus", cascade="all, delete-orphan"))

    __table_args__ = (
        # UniqueConstraint(gene_id, taxon_id, name=__tablename__ + '_c1'),
        UniqueConstraint(taxon_id, left, right, name=__tablename__ + '_c1'),
    )


class Loci(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}locis"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))


class LociLocus(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}loci_locus"

    loci_id = Column(BigInteger, ForeignKey(Loci.id), nullable=False, primary_key=True)
    loci = relationship(Loci, backref=backref("loci", cascade="all, delete-orphan"))
    locus_id = Column(BigInteger, ForeignKey(Locus.id), nullable=False, primary_key=True)
    locus = relationship(Locus, backref=backref("loci", cascade="all, delete-orphan"))


class Specimen(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}specimens"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['specimen'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id, ondelete='CASCADE'), primary_key=True)


class Sequence(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}sequences"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['sequence'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id, ondelete='CASCADE'), primary_key=True)
    specimen_id = Column(BigInteger, ForeignKey(Specimen.id), nullable=True, primary_key=False)
    specimen = relationship(Specimen, backref=backref("sequences", cascade="all, delete-orphan"),
                            foreign_keys=[specimen_id])


class MultipleSequenceAlignment(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}msas"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['multiple-sequence-alignment'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id, ondelete='CASCADE'), primary_key=True)


class PhylogeneticTree(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}phylo_trees"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['phylogenetic-tree'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id, ondelete='CASCADE'), primary_key=True)


class SequenceSimilarity(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}seq_sims"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['sequence-similarity'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id, ondelete='CASCADE'), primary_key=True)


class DiscriminantMatrix(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}disc_matrix"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['discriminant-matrix'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id, ondelete='CASCADE'), primary_key=True)


class Supermatrix(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}supermatrix"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['supermatrix'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id, ondelete='CASCADE'), primary_key=True)
