import pandas as pd
import geopandas as gpd
from geo.Geoserver import Geoserver
from flask import Blueprint, request, g
from flask.views import MethodView
from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject, Issue, IType, register_api
from biobarcoding.db_models.geographics import GeographicRegion, GeographicLayer, Regions
import json
from biobarcoding.geo import geoserver_session
from biobarcoding import postgis_engine


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

def get_content(session, Feature, id=None):
    content = None
    count = 0
    try:
        content = session.query(Feature)
        if id:
            content = content.filter(Feature.id == id).first()
        else:
            content = content.order_by(Feature.id).all()
        issues, status = [Issue(IType.INFO, f'READ "{Feature.__name__}": The "{Feature.__name__}" were successfully read')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'READ "{Feature.__name__}": The "{Feature.__name__}" could not be read.')], 500
    if content == None:
        issues, status = [Issue(IType.ERROR, f'no data available')], 200
        content = ""
    return issues, content, count, status

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

def get_json_from_schema(entity, t):
    entity_schema = getattr(entity, "Schema")()
    t_json = entity_schema.dumps(t)
    return json.loads(t_json)



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
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        issues, geographicregion, count, status = get_content(db, GeographicRegion,region_id)
        if status == 200 and geographicregion != "":
            if region_id is None:
                lines = False
                issues, regions, count, status = get_content(pg, Regions,region_id)
            else:
                lines = True
                issues, regions, count, status = get_content(pg, Regions,geographicregion.geo_id)

            geodf = pd.read_json(response_to_dataframe(geographicregion),lines = lines)
            df = pd.read_json(response_to_dataframe(regions), lines = lines)
            geodf = geodf.set_index(geodf["geo_id"]).drop(columns=["geo_id","uuid","id"])
            df = df.set_index(df["id"]).drop(columns=["id","uuid"])
            content = pd.concat([geodf,df], axis = 1, join="inner")
        else:
            content  = None
        return ResponseObject(issues= issues,status = status, content= content, content_type= "application/json").get_response()


    @bcs_session()
    def post(self):
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        r = ResponseObject()
        t = request.json
        GeographicRegion_schema = getattr(GeographicRegion, "Schema")()
        Regions_schema = getattr(Regions, "Schema")()
        regions_data = get_json_from_schema(Regions, t)
        geographicregion_data = get_json_from_schema(GeographicRegion, t)
        regions = Regions_schema.load(regions_data, instance=Regions())
        pg.add(regions)
        pg.commit()
        geographicregion = GeographicRegion_schema.load(geographicregion_data, instance = GeographicRegion())
        geographicregion.uuid = regions.uuid
        geographicregion.geo_id = regions.id
        geographicregion.usr = g.bcs_session.identity.id
        db.add(geographicregion)
        return r.get_response()


    @bcs_session()
    def delete(self,region_id = None):
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        issues, geographicregion, count, status = get_content(db, GeographicRegion, region_id)
        if status == 200 and geographicregion != "":
            issues, region, count, status = get_content(pg, Regions,geographicregion.geo_id)
            db.delete(geographicregion)
            pg.delete(region)
        return ResponseObject(issues= issues,status = status).get_response()


    @bcs_session()
    def put(self, region_id = None):
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        t = request.json
        GeographicRegion_schema = getattr(GeographicRegion,"Schema")()
        Regions_schema = getattr(Regions, "Schema")()
        issues, geographicregion, count, status = get_content(db, GeographicRegion, region_id)
        if status == 200 and geographicregion != "":
            issues, region, count, status = get_content(pg, Regions, geographicregion.geo_id)
        geographicregion = GeographicRegion_schema.load(get_json_from_schema(GeographicRegion, t), instance = geographicregion)
        regions = Regions_schema.load(get_json_from_schema(Regions, t), instance = region)
        db.add(geographicregion)
        pg.add(regions)
        return ResponseObject(issues= issues,status = status).get_response()


register_api(bp_geo, RegionsAPI, "geo/regions", f"{bcs_api_base}/geo/regions/", pk="region_id")

class layerAPI():
    '''
    geo_session... inicializar en main o inicializar cada vez que se haga una llamada restful?)
    GET: - the list of layers
         - publish a result of filtering a layer (with geoserver_session.publish_featurestore_sqlview(sql_string))

    POST: create and store (in postgis) and publish a new layer in which format am I going to retyreve that vector layer?
    PUT: changing attributes
    DELETE:
    '''
    @bcs_session()
    def get(self, id = None):
        db = g.bcs_session.db_session
        issues, layers, count, status = get_content(db, GeographicLayer, id)
        return ResponseObject(issues= issues,status = status, content= layers)

    @bcs_session()
    def post(self):
        '''
        json with
        name:''
        geopandas:''/filter:''
        @return:
        '''
        db = g.bcs_session.db_session
        t = request.json
        GeographicLayer_Schema= getattr(GeographicLayer,"Schema")()
        geographiclayer = GeographicLayer_Schema.load(t, instance=GeographicLayer())
        db.add(geographiclayer)
        df = gpd.read(t["geojson"])
        try:
            df.to_postgis(t["name"],postgis_engine, if_exists="fail") # names have to be uniques or will be stores by id
        except ValueError:
            issues, status = [Issue(IType.INFO,f"layer named as {t['name']} already exists")], 400

        else:
            geoserver_session.publish_featurestore(workspace='ngd', store_name='geo_data', pg_table=t['name'])
            issues, status = [Issue(IType.INFO,f"layer named as {t['name']} published")], 200
        return ResponseObject(issues=issues, status=status, content=None)

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
    pass

def publish_tmp_layer(gdf, workspace, data_storage_name):
    geo = Geoserver('http://127.0.0.1:8080/geoserver', username='admin', password='geoserver')
    pass


def mask_raster_from_layer(raster, vector):
    geo = Geoserver('http://127.0.0.1:8080/geoserver', username='admin', password='geoserver')
    raster_layer = geo.get_layer(layer_name=raster)
    vector_layer = geo.get_layer(layer_name=vector)
    print(vector_layer)
    geo.create_featurestore(store_name='geo_data', workspace='ngd', db='ngd_geoserver', host='localhost', port=5435,
                            pg_user='postgres', pg_password='postgres')
    pass