# GEOGRAPHIC data
import uuid

from geoalchemy2 import Geometry
from sqlalchemy import Column, Integer, String, BigInteger, ForeignKey, JSON, Boolean, Index
from sqlalchemy.dialects.postgresql import JSONB

from . import ORMBase, ORMBaseGeo, GUID
from .core import Dataset, data_object_type_id

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


class GeographicLayer(Dataset):
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['geolayer'],
    }
    __tablename__ = f"{prefix}layers"

    id = Column(BigInteger, ForeignKey(Dataset.id), primary_key=True)

    # NOTE: identity_id, name and attributes, in "dataset"
    # identity_id = Column(Integer)
    # name = Column(String(80))
    # attributes = Column(JSON)  # Metadata about the layer: categories, tags, ...
    wks = Column(String(80))  # Workspace
    geoserver_name = Column(String(80))
    wms_url = Column(String(255))
    properties = Column(JSONB)  # Store information about fields in a vectorial layer. Name, style, range, ...
    published = Column(Boolean, default=False)
    in_postgis = Column(Boolean, default=False)
    is_deleted = Column(Boolean, default=False)
    layer_type = Column(String(80))  # CHOICES?


Index("index_GeographicLayer_on_Properties_gin",
      GeographicLayer.properties,
      postgresql_using="gin")


# NOTE!!! This table is a PostGIS table (not in App DB; notice the parent class)
# It is a special feature layer used to contain the geometry of regions
class Regions(ORMBaseGeo):
    __tablename__ = 'Regions'

    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    id = Column(Integer, primary_key=True, autoincrement=True)
    geometry = Column(Geometry(geometry_type='MULTIPOLYGON', srid=4326))
