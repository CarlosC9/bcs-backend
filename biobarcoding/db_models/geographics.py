# GEOGRAPHIC data
from sqlalchemy import Column, Integer, String

from biobarcoding.db_models import ORMBase, GUID

prefix = "geo_"

class GeographicRegion(ORMBase):
    __tablename__ = f"{prefix}regions"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class GeographicLayer(ORMBase):
    __tablename__ = f"{prefix}layers"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
