import pandas as pd
from sqlalchemy import create_engine
import geopandas as gpd
from geo.Geoserver import Geoserver
from flask import Blueprint, request, send_file, g, make_response, jsonify
from flask.views import MethodView
from biobarcoding.authentication import bcs_session
from biobarcoding.common import generate_json
from biobarcoding.rest import bcs_api_base, bcs_gui_base, ResponseObject, Issue, IType, register_api
from biobarcoding.db_models import DBSession, DBSessionGeo
from biobarcoding.db_models.geographics import GeographicRegion, GeographicLayer, Regions
from shapely.geometry import shape
from shapely.geometry.multipolygon import MultiPolygon
import json
from geoalchemy2.shape import to_shape
from sqlalchemy.orm import sessionmaker, scoped_session
import os

"""
Support operations regarding Geographical layers

* Store raster layer in both:
  - Geoserver. For display
  - PostGIS. For calculations?

* How to translate temporary-user layers to unique PostGIS layers? (both feature and raster)
  - Postgres Schemas??
  - Or translate a GUI query to a PostGIS query?

New ENTITY "geolayer"
* tags/attributes.

Four areas:
* system wide layers
* system wide regions
* user layers
* user regions

* CRUD feature layer
  - list
  - read (WMS, WFS)
  - import. (name, ¿projection?, workspace, user, format, [sld], store (Geoserver or PostGIS), <file binary>)
  - export
* CRUD import export coverage layer. Includes PDA layers
* CRUD SLDs
  - select SLD for a layer. Default SLD for layer
  - generate SLD for common colormaps
* CRUD regions. A region is a multipolygon which can be used to filter
  - A feature table (PostGIS) for system REGIONS
  - A feature table per user for regions
  - How a region is created?
    - Direct
    - Direct selection of one or more existing polygons
    - Output of a GEOS function, which can include features and regions
      - A special type of function whose result would be a geopandas with no value columns?
        - The geometry would be the ¿OR/AND? combination of all rows
* Setup NGD workspaces: NGD-SYS, and NGD-USR
* Setup NGD layer groups: a tree, then assign layers into one (or more?) of them
  - PDA
  - Biota
  - Grafcan
  -
* Output raster layer. Apply SLD automatically, store in users workspace
* Output feature layer. Apply SLD automatically, store in users workspace
"""

bp_geo = Blueprint('geo', __name__)


# def inizialize_regions_table():
#     session = DBSessionGeo()
#     dummy_geojson = {"features": [
#         {
#             "type": "Feature",
#             "properties": {},
#             "geometry": {
#                 "type": "Polygon",
#                 "coordinates": [
#                     [
#                         [
#                             -15.503768920898436,
#                             28.232414812316087
#                         ],
#                         [
#                             -15.500335693359375,
#                             28.043500801345925
#                         ],
#                         [
#                             -15.303955078125,
#                             28.00531432138343
#                         ],
#                         [
#                             -15.292968749999998,
#                             28.168875180063345
#                         ],
#                         [
#                             -15.503768920898436,
#                             28.232414812316087
#                         ]
#                     ]
#                 ]
#             }
#         }
#     ]
#     }
#     geom = MultiPolygon([shape(i['geometry']) for i in dummy_geojson["features"]])
#     region = Regions(name='dummy_region', geometry=f'SRID=4326; {geom.wkt}')
#     session.add(region)
#     session.commit()



class RegionsAPI(MethodView):
    '''
    '''
    @bcs_session()
    def get(self, region_id = None):
        '''
        ?? dede la GUI recibo el polígono marcado por el usuario
        @return:
        # dar lista de regiones para visualizar
        # entregar región a PDA
        '''
        # tengo que combinar con la otra tabla
        db = g.bcs_session.db_session
        geo_session = DBSessionGeo()
        rgeo = ResponseObject()
        rr = ResponseObject()
        if region_id is None:
            lines = False
            rgeo.content = db.query(GeographicRegion).order_by(GeographicRegion.name).all()
            res_geo = rgeo.get_response()
            if rgeo.status == 200:
                rr.content = geo_session.query(Regions).order_by(Regions.name).all()
        else:
            lines = True
            rgeo.content = db.query(GeographicRegion).filter(GeographicRegion.id == region_id).first()
            res_geo = rgeo.get_response()
            if rgeo.status == 200:
                rr.content = geo_session.query(Regions).filter(Regions.id == rgeo.content.geo_id).first()
        if rr.status == 200:
            geodf = pd.read_json(response_to_dataframe(res_geo),lines = lines)
            df = pd.read_json(response_to_dataframe(rr.get_response()), lines = lines)
            geodf = geodf.set_index(geodf["geo_id"]).drop(columns=["geo_id"])
            df = df.set_index(df["id"]).drop(columns=["id"])
            df = pd.concat([geodf,df], axis = 1, join="inner")
            return make_response(jsonify(df.to_dict("index"))), 200
        else:
            return ResponseObject(status=200)


    # @bcs_session()
    def post(self):
        '''

        @param poligon: geoJSON format
        @param name: name of region
        @return: response
        '''
        # Añadir región a tabla en postgis
        # Añadir region a GeographicLayer
        t = request.json
        geojson = t.get('geojson')
        geom = MultiPolygon([shape(i['geometry']) for i in geojson["features"]])
        session = DBSessionGeo()
        region = Regions(name='dummy_region', geometry=f'SRID=4326; {geom.wkt}')
        session.add(region)
        session.commit()
        # actualizar
        geograficregion = GeographicRegion()
        session = DBSession()
        geograficregion.name = t.get('name')
        geograficregion.usr = t.get('usr')
        geograficregion.tags = t.get('tags')
        geograficregion.geo_id = region.id
        session.add(geograficregion)
        session.commit()
        response_object = {
            'status': 'success',
            'message': f"region with ID {geograficregion.id} addaed to layer Regions"
        }
        DBSessionGeo.remove()
        return make_response(jsonify(response_object)), 200

    # @bcs_session
    def delete(self,region_id = None):
        session = DBSession()
        geosession = DBSessionGeo()
        r = ResponseObject()
        geographicregion = session.query(GeographicRegion).filter(GeographicRegion.id == region_id).first()
        region = geosession.query(Regions).filter(Regions.id == geographicregion.geo_id).first()
        session.delete(geographicregion)
        geosession.delete(region)
        session.commit()
        geosession.commit()
        return r.get_response()


    # @bcs_session
    def put(self, region_id = None):
        session = DBSession()
        GEOsession = DBSessionGeo()
        t = request.json
        GeographicRegion_schema = getattr(GeographicRegion,"Schema")()
        Regions_schema = getattr(Regions, "Schema")()
        geographicregion = session.query(GeographicRegion).filter(GeographicRegion.id == region_id).first()
        regions = GEOsession.query(Regions).filter(Regions.id == geographicregion.id).first()
        geographicregion = GeographicRegion_schema.load(self._prepare(GeographicRegion,t),instance = geographicregion)
        regions = Regions_schema.load(self._prepare(Regions,t),instance = regions)
        session.add(geographicregion)
        GEOsession.add(regions)
        GEOsession.commit()
        session.commit()
        return jsonify(message='Successfuly updated'), 200


    def _prepare(self,entity,t):
        entity_schema = getattr(entity,"Schema")()
        t_json = entity_schema.dumps(t)
        return json.loads(t_json)

register_api(bp_geo, RegionsAPI, "geo/regions", f"{bcs_api_base}/geo/regions/", pk="region_id")


def response_to_dataframe(response):
    data = response.get_data()
    my_json = data.decode('utf8').replace("'", '"')
    data = json.loads(my_json)
    content = data.get("content")
    data = json.dumps(content)
    return data
    # lines = True para un item
    # lines = False para varios



class layerAPI():
    @bcs_session()
    def get(self, id, filter):
        pass

    @bcs_session()
    def post(self, id):
        pass

    @bcs_session()
    def put(self, id):
        pass

    @bcs_session()
    def delete(self, id):
        '''
        delete layer from PostGIS
        ¿and geoserver in case of raster layer?
        @param id:
        @return:
        '''
        pass


class stylesAPI():

    @bcs_session()
    def get(self, id, layer, style_name):
        '''
        apply an style in a layer and publish
        @param id:
        @return:
        '''
        # geo.publish_style(layer_name=layer, style_name='sld_file_name', workspace=style_name,
        #                 sld_version='1.0.0')# version?
        pass

    @bcs_session()
    def put(self, rampa):
        '''
        create a new stylein geoserver
        @param rampa:
        @return:
        '''
        # create new style
        # geo.upload_style(path=r'path\to\sld\file.sld', workspace='demo')
        pass

def new_layer_data_filter(layer, filter):
    geo = Geoserver('http://127.0.0.1:8080/geoserver', username='admin', password='geoserver')
    sql = 'SELECT geometry as geom, "DENOMTAX" FROM plantas WHERE "RAREZALOCA" > 1'
    geo.publish_featurestore_sqlview(store_name='geo_data', name='tmp_view3', sql=sql, key_column="'IDCELDA'",
                                     workspace='ngd')

def publish_tmp_layer(gdf, workspace, data_storage_name):
    geo = Geoserver('http://127.0.0.1:8080/geoserver', username='admin', password='geoserver')


def mask_raster_from_layer(raster, vector):
    geo = Geoserver('http://127.0.0.1:8080/geoserver', username='admin', password='geoserver')
    # raster_layer = geo.get_layer(layer_name=raster)
    vector_layer = geo.get_layer(layer_name=vector)
    print(vector_layer)
    pass
