import pandas as pd
from geo.Geoserver import Geoserver
from flask import Blueprint, request, send_file, g, make_response, jsonify
from flask.views import MethodView
from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, bcs_gui_base, ResponseObject, Issue, IType, register_api
from biobarcoding.db_models import DBSession, DBSessionGeo
from biobarcoding.db_models.geographics import GeographicRegion, GeographicLayer, Regions
import json


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

class RegionsAPI(MethodView):
    '''
    '''
    @bcs_session(read_only=True)
    def get(self, region_id = None):
        '''
        ?? dede la GUI recibo el polígono marcado por el usuario
        @return:
        # dar lista de regiones para visualizar
        # entregar región a PDA
        '''
        # tengo que combinar con la otra tabla
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        if region_id is None:
            lines = False
            geographicregion = db.query(GeographicRegion).order_by(GeographicRegion.name).all()
            regions = pg.query(Regions).order_by(Regions.name).all()
        else:
            lines = True
            geographicregion = db.query(GeographicRegion).filter(GeographicRegion.id == region_id).first()
            regions = pg.query(Regions).filter(Regions.id == geographicregion.geo_id).first()

        geodf = pd.read_json(response_to_dataframe(geographicregion),lines = lines)
        df = pd.read_json(response_to_dataframe(regions), lines = lines)
        geodf = geodf.set_index(geodf["geo_id"]).drop(columns=["geo_id"])
        df = df.set_index(df["id"]).drop(columns=["id"])
        df = pd.concat([geodf,df], axis = 1, join="inner")
        return make_response(jsonify(df.to_dict("index"))), 200



    @bcs_session()
    def post(self):
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        r = ResponseObject()
        t = request.json
        GeographicRegion_schema = getattr(GeographicRegion, "Schema")()
        Regions_schema = getattr(Regions, "Schema")()
        regions_data = self._get_json_from_schema(Regions, t)
        geographicregion_data = self._get_json_from_schema(GeographicRegion, t)
        regions = Regions_schema.load(regions_data, instance=Regions())
        pg.add(regions)
        pg.commit()
        # TODO se debería forzar a error si no se aporta geometría?
        geographicregion_data["geo_id"] = regions.id
        geographicregion = GeographicRegion_schema.load(geographicregion_data, instance = GeographicRegion())
        db.add(geographicregion)
        return r.get_response()


    @bcs_session()
    def delete(self,region_id = None):
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        r = ResponseObject()
        geographicregion = db.query(GeographicRegion).filter(GeographicRegion.id == region_id).first()
        region = pg.query(Regions).filter(Regions.id == geographicregion.geo_id).first()
        db.delete(geographicregion)
        pg.delete(region)
        return r.get_response()


    @bcs_session()
    def put(self, region_id = None):
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        r = ResponseObject()
        t = request.json
        GeographicRegion_schema = getattr(GeographicRegion,"Schema")()
        Regions_schema = getattr(Regions, "Schema")()
        geographicregion = db.query(GeographicRegion).filter(GeographicRegion.id == region_id).first()
        regions = pg.query(Regions).filter(Regions.id == geographicregion.id).first()
        geographicregion = GeographicRegion_schema.load(self._get_json_from_schema(GeographicRegion, t), instance = geographicregion)
        regions = Regions_schema.load(self._get_json_from_schema(Regions, t), instance = regions)
        db.add(geographicregion)
        pg.add(regions)
        return r.get_response()


    def _get_json_from_schema(self,entity,t):
        entity_schema = getattr(entity,"Schema")()
        t_json = entity_schema.dumps(t)
        return json.loads(t_json)

register_api(bp_geo, RegionsAPI, "geo/regions", f"{bcs_api_base}/geo/regions/", pk="region_id")


def response_to_dataframe(item):
    r = ResponseObject()
    r.content = item
    res = r.get_response()
    data = res.get_data()
    my_json = data.decode('utf8').replace("'", '"')
    data = json.loads(my_json)
    content = data.get("content")
    data = json.dumps(content)
    return data



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
