# BIOINFORMATIC data
from sqlalchemy import Column, Integer, ForeignKey, UnicodeText, Unicode, String, BigInteger, UniqueConstraint
from sqlalchemy.orm import relationship

from biobarcoding.db_models import ORMBase, GUID, ObjectType
import uuid

prefix = "bo_"

bio_object_type_id = {
    "sequence": 1,
    "multiple-sequence-alignment": 2,
    "phylogenetic-tree": 3,
    "geographic-layer": 4,
    "sequence-similarity": 5,  # BLAST
    "dataframe": 6
}


class BioinformaticObject(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}bos"
    id = Column(BigInteger, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    bo_type_id = Column(Integer, ForeignKey(ObjectType.id))
    chado_id = Column(BigInteger, nullable=False)  # Foreign key (not enforceable by the DB)
    chado_table = Column(String(80))
    name = Column(String(80))
    # content = Column(UnicodeText)

    __table_args__ = (
        UniqueConstraint(chado_id, chado_table, name=__tablename__+'_c1'),
    )
    __mapper_args__ = {
        'polymorphic_identity': 'bioinformatic_obj',
        'polymorphic_on': bo_type_id
    }


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


class Sequence(BioinformaticObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}sequences"
    __mapper_args__ = {
        'polymorphic_identity': bio_object_type_id['sequence'],
    }
    id = Column(BigInteger, ForeignKey(BioinformaticObject.id), primary_key=True)
    specimen_id = Column(Integer, ForeignKey(Specimen.id), nullable=True, primary_key=False)
    specimen = relationship(Specimen)


class MultipleSequenceAlignment(BioinformaticObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}msas"
    __mapper_args__ = {
        'polymorphic_identity': bio_object_type_id['multiple-sequence-alignment'],
    }
    id = Column(BigInteger, ForeignKey(BioinformaticObject.id), primary_key=True)


class PhylogeneticTree(BioinformaticObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}phylo_trees"
    __mapper_args__ = {
        'polymorphic_identity': bio_object_type_id['phylogenetic-tree'],
    }
    id = Column(BigInteger, ForeignKey(BioinformaticObject.id), primary_key=True)


class SequenceSimilarity(BioinformaticObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}seq_sims"
    __mapper_args__ = {
        'polymorphic_identity': bio_object_type_id['sequence-similarity'],
    }
    id = Column(BigInteger, ForeignKey(BioinformaticObject.id), primary_key=True)
