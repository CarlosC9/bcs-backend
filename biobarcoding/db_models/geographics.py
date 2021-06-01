# GEOGRAPHIC data
from marshmallow_sqlalchemy import SQLAlchemySchema, fields
from sqlalchemy import Column, Integer, String, BigInteger, ForeignKey, JSON, Boolean

from biobarcoding.db_models import ORMBase, ORMBaseGeo, GUID
from biobarcoding.db_models.bioinformatics import BioinformaticObject, bio_object_type_id
from geoalchemy2 import Geometry
import uuid

prefix = "geo_"


class GeographicRegion(ORMBase):
    __tablename__ = f"{prefix}regions"

    uuid = Column(GUID)
    id = Column(Integer, primary_key=True, autoincrement=True)
    geo_id = Column(Integer, unique=True)
    name = Column(String(80))
    usr = Column(Integer)
    attributes = Column(JSON)


class GeographicLayer(ORMBase):
    __tablename__ = f"{prefix}layers"

    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    id = Column(Integer, primary_key=True, autoincrement=True)
    usr = Column(Integer)
    name = Column(String(80))
    wks = Column(String(80))
    attributes = Column(JSON)
    published = Column(Boolean,default= False)
    in_postgis = Column(Boolean,default= False)

class Regions(ORMBaseGeo):
     __tablename__ = 'Regions'

     uuid = Column(GUID, unique=True, default=uuid.uuid4)
     id = Column(Integer, primary_key=True, autoincrement=True)
     geometry = Column(Geometry(geometry_type='MULTIPOLYGON', srid=4326))




