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
    identity_id = Column(Integer)
    attributes = Column(JSON)
    style = Column(JSON)


class GeographicLayer(ORMBase):
    __tablename__ = f"{prefix}layers"

    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    id = Column(Integer, primary_key=True, autoincrement=True)
    identity_id = Column(Integer)
    name = Column(String(80))
    wks = Column(String(80))  # Workspace
    attributes = Column(JSON)  # Metadata about the layer: categories, tags, ...
    geoserver_name = Column(String(80))
    wms_url = Column(String(255))
    properties = Column(JSON)  # Store information about fields in a vectorial layer. Name, style, range, ...
    published = Column(Boolean, default=False)
    in_postgis = Column(Boolean, default=False)
    is_deleted = Column(Boolean, default=False)
    layer_type = Column(String(80))  # CHOICES?


# This table is a PostGIS table (not in BCS; notice the parent class)
# It is a special feature layer used to contain the geometry of regions
class Regions(ORMBaseGeo):
     __tablename__ = 'Regions'

     uuid = Column(GUID, unique=True, default=uuid.uuid4)
     id = Column(Integer, primary_key=True, autoincrement=True)
     geometry = Column(Geometry(geometry_type='MULTIPOLYGON', srid=4326))




