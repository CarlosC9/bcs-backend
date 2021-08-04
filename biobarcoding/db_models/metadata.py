import datetime

from sqlalchemy import Column, Integer, String, Boolean, DateTime

from biobarcoding.db_models import ORMBase, GUID

prefix = "md_"


class Ontology(ORMBase):
    __tablename__ = f"{prefix}ontologies"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class Publication(ORMBase):
    __tablename__ = f"{prefix}publications"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class Taxonomy(ORMBase):
    __tablename__ = f"{prefix}taxonomies"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class Taxon(ORMBase):
    __tablename__ = f"{prefix}taxa"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class SpeciesNameToCanonical(ORMBase):
    __tablename__ = f"{prefix}species_names"

    name = Column(String(500), primary_key=True)
    rank = Column(String(20))
    synonym = Column(Boolean)
    canonical_name = Column(String(200))
    scientific_name = Column(String(200))
    is_canonical = Column(Boolean, default=False)
    creation_time = Column(DateTime, default=datetime.datetime.utcnow())
