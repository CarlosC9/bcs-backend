# GEOGRAPHIC data
from sqlalchemy import Column, Integer, String, BigInteger, ForeignKey

from biobarcoding.db_models import ORMBase, GUID
from biobarcoding.db_models.bioinformatics import BioinformaticObject, bio_object_type_id

prefix = "geo_"


class GeographicRegion(ORMBase):
    __tablename__ = f"{prefix}regions"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class GeographicLayer(BioinformaticObject):
    __tablename__ = f"{prefix}layers"
    __mapper_args__ = {
        'polymorphic_identity': bio_object_type_id['geographic-layer'],
    }
    id = Column(BigInteger, ForeignKey(BioinformaticObject.id), primary_key=True)
