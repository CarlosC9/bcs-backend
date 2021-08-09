# GEOGRAPHIC data
import uuid

from geoalchemy2 import Geometry
from sqlalchemy import Column, Integer, String, JSON, Boolean
from sqlalchemy.dialects.postgresql import JSONB

from . import ORMBase, ORMBaseGeo, GUID

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
    # TODO To migrate existing "geo_layers" table:
    #  alter table geo_layers alter column properties type jsonb using properties::jsonb;
    #  alter table geo_layers alter column attributes type jsonb using attributes::jsonb;
    attributes = Column(JSONB)  # Metadata about the layer: categories, tags, ...
    properties = Column(JSONB)  # Store information about fields of a vectorial layer: field name, type, style, range...
    geoserver_name = Column(String(80))
    wms_url = Column(String(255))
    published = Column(Boolean, default=False)
    in_postgis = Column(Boolean, default=False)
    is_deleted = Column(Boolean, default=False)
    layer_type = Column(String(80))  # CHOICES?


# NOTE!!! This table is a PostGIS table (not in BCS; notice the parent class)
# It is a special feature layer used to contain the geometry of regions
class Regions(ORMBaseGeo):
    __tablename__ = 'Regions'

    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    id = Column(Integer, primary_key=True, autoincrement=True)
    geometry = Column(Geometry(geometry_type='MULTIPOLYGON', srid=4326))
