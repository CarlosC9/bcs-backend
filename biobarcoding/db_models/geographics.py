# GEOGRAPHIC data
from marshmallow_sqlalchemy import SQLAlchemySchema, fields
from sqlalchemy import Column, Integer, String, BigInteger, ForeignKey, JSON, MetaData

from biobarcoding.db_models import ORMBase, ORMBaseGeo, GUID
from biobarcoding.db_models.bioinformatics import BioinformaticObject, bio_object_type_id
from geoalchemy2 import Geometry
from geoalchemy2.shape import to_shape

prefix = "geo_"


class GeographicRegion(ORMBase):
    __tablename__ = f"{prefix}regions"

    id = Column(Integer, primary_key=True, autoincrement=True)
    geo_id = Column(Integer, unique=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    usr = Column(String(80))
    tags = Column(JSON)


class GeographicLayer(BioinformaticObject):
    __tablename__ = f"{prefix}layers"
    __mapper_args__ = {
        'polymorphic_identity': bio_object_type_id['geographic-layer'],
    }
    id = Column(BigInteger, ForeignKey(BioinformaticObject.id), primary_key=True)


class Regions(ORMBaseGeo):
     __tablename__ = 'Regions'
     id = Column(Integer, primary_key=True)
     name = Column(String)
     geometry = Column(Geometry(geometry_type='MULTIPOLYGON', srid=4326))




