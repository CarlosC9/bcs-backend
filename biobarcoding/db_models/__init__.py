import json
import uuid
from copy import deepcopy

import geopandas
from geoalchemy2.elements import WKBElement
from geoalchemy2.shape import to_shape
from geoalchemy2.types import Geometry as Multipolygon
from marshmallow.fields import String as SchemaString
from marshmallow_sqlalchemy import ModelConversionError, ModelSchema
from marshmallow_sqlalchemy import ModelConverter as BaseModelConverter
from shapely.geometry import shape
from shapely.geometry.multipolygon import MultiPolygon
from sqlalchemy import event, TypeDecorator, CHAR, Column, Integer, String
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker, class_mapper, ColumnProperty, RelationshipProperty, mapper
from sqlalchemy_continuum import make_versioned

# make_versioned(user_cls=None, options={'native_versioning': True})
make_versioned(user_cls=None)

DBSession = scoped_session(sessionmaker())
DBSessionChado = scoped_session(sessionmaker())
DBSessionGeo = scoped_session(sessionmaker())


class GUID(TypeDecorator):
    """
    Platform-independent GUID type.

    Uses Postgresql's UUID type, otherwise uses
    CHAR(32), storing as stringified hex values.

    """
    impl = CHAR

    def load_dialect_impl(self, dialect):
        if dialect.name == 'postgresql':
            return dialect.type_descriptor(UUID())
        else:
            return dialect.type_descriptor(CHAR(32))

    def process_bind_param(self, value, dialect):
        if value is None:
            return value
        elif dialect.name == 'postgresql':
            return str(value)
        else:
            if not isinstance(value, uuid.UUID):
                return "%.32x" % uuid.UUID(value).int
            else:
                # hexstring
                return "%.32x" % value.int

    def process_result_value(self, value, dialect):
        # if dialect.name == 'postgresql':
        #     return uuid.UUID(value)
        # else:
        #     return value
        if value is None:
            return value
        else:
            return uuid.UUID(value)


class BaseMixin(object):
    # query = DBSession.query_property()
    # @declared_attr
    # def __tablename__(cls):
    # return cls.__name__.lower()

    # def __new__(cls, *args, **kwargs):
    #     obj = super().__new__(cls)
    #     if '_sa_instance_state' not in obj.__dict__:
    #         obj._sa_instance_state = InstanceState(obj, obj._sa_class_manager)
    #     return obj

    def __deepcopy__(self, memo):
        cls = self.__class__
        print(str(cls))
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k in class_mapper(cls).iterate_properties:
            if isinstance(k, ColumnProperty):
                name = k.columns[0].name
                getattr(self, name, None)
            elif isinstance(k, RelationshipProperty):
                # if k.back_populates:
                name = k.strategy.key
                getattr(self, name, None)

        for k, v in self.__dict__.items():
            deepcopy(v, memo)
            # setattr(result, k, deepcopy(v, memo))
        return result


ORMBase = declarative_base(cls=BaseMixin)
ORMBaseChado = declarative_base(cls=BaseMixin)
ORMBaseGeo = declarative_base(cls=BaseMixin)


class ObjectType(ORMBase):  # CODES
    """ Sequence, Alignment, Phylogenetic Tree, ... but also "Functions of the system" (like import, export, etc.) """
    __tablename__ = "object_types"

    id = Column(Integer, primary_key=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80), nullable=False)


class GeoString(SchemaString):
    """A string field.

    :param kwargs: The same keyword arguments that :class:`Field` receives.
    """

    #: Default error messages.
    default_error_messages = {
        "invalid": "Not a valid string.",
        "invalid_utf8": "Not a valid utf-8 string.",
    }

    def _serialize(self, value, attr, obj, **kwargs):
        if value is None:
            return None
        if not isinstance(value, (dict, WKBElement)):
            raise self.make_error("invalid")
        if isinstance(value, WKBElement):
            multi = to_shape(value)
            polygons = list(multi)
            geojson = geopandas.GeoSeries(polygons).to_json()
            geodict = json.loads(geojson)
            # parse geodict to get rid og bbox attribute.
            # I assume that it will be a bbox attribute for feature and one bbox for the Feature collection
            # I assume that Featureccollection is a list of Polygones
            geodict.pop('bbox', None)
            for geo in geodict['features']:
                geo.pop('bbox', None)
            return geodict
        if isinstance(value, dict):
            return value

    def _deserialize(self, value, attr, data, **kwargs):
        if not isinstance(value, (dict, bytes)):
            raise self.make_error("invalid")
        else:
            geom = MultiPolygon([shape(i['geometry']) for i in value["features"]])
            return f'SRID=4326; {geom.wkt}'


def setup_schema(Base, session):
    from biobarcoding.db_models.geographics import Regions

    # Create a function which incorporates the Base and session information

    def setup_schema_fn():
        # https://stackoverflow.com/questions/58299923/serialize-geometry-using-flask-marshmallow
        # https://github.com/marshmallow-code/marshmallow-sqlalchemy/issues/55

        class ModelConverter(BaseModelConverter):
            SQLA_TYPE_MAPPING = {
                **BaseModelConverter.SQLA_TYPE_MAPPING,
                **{Multipolygon: GeoString},
            }

        class GeoSchema(ModelSchema):
            class Meta:
                model = Regions
                model_converter = ModelConverter

        for class_ in Base._decl_class_registry.values():
            if hasattr(class_, "__tablename__"):
                if class_.__name__.endswith("Schema"):
                    raise ModelConversionError(
                        "For safety, setup_schema can not be used when a"
                        "Model class ends with 'Schema'"
                    )

                class Meta(object):
                    include_fk = True
                    if hasattr(class_, "_ts_vector"):
                        exclude = ("_ts_vector",)
                    model = class_
                    sqla_session = session
                    generate_polymorphic_schemas = True

                schema_class_name = "%sSchema" % class_.__name__

                if class_.__name__ == "Regions":
                    schema_class = type(schema_class_name, (GeoSchema,), {"Meta": Meta})
                else:
                    schema_class = type(schema_class_name, (ModelSchema,), {"Meta": Meta})

                setattr(class_, "Schema", schema_class)

    return setup_schema_fn


event.listen(mapper, "after_configured", setup_schema(ORMBase, DBSession))
event.listen(mapper, "after_configured", setup_schema(ORMBaseChado, DBSessionChado))
event.listen(mapper, "after_configured", setup_schema(ORMBaseGeo, DBSessionGeo))
