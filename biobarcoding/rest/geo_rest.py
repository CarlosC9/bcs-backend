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
import pathlib
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT


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

class layerAPI(MethodView):
    '''
    geo_session... inicializar en main o inicializar cada vez que se haga una llamada restful?)
        GET: - the list of layers info from Geographiclayer table on bcs
             - publish a result of filtering a layer (with geoserver_session.publish_featurestore_sqlview(sql_string))

    POST:   Three post cases (at least):
            1. post all new vector layer from geoJSON format (PDA case)
            2. Import a Raster o Vector layer from a file (via command)
            3. Save a temporal layer created by geoserver_session.publish_featurestore_sqlview(sql_string)
            create and store in geographichlayer table as item
            then in postgis as table (format)
            and then publish  in geoserver as new layer

    PUT: changing name or attribues (never Geometry)

    DELETE: delete completely from GeographicLayer table (bcs) postgis and geoserver

    '''
    @bcs_session()
    def get(self, filter = None):
        '''
        - get the list of layers in postgis
        - create a temporal layer from a sql filter
        @param filter: SQL filter in case of temporal layer
        @return:
        '''
        if filter == None:
            db = g.bcs_session.db_session
            issues, layers, count, status = get_content(db, GeographicLayer)
            return ResponseObject(issues=issues, status=status, content=layers).get_response()
        else:
            # el tema de la key_columns...???
            r = geoserver_session.publish_featurestore_sqlview(store_name='geo_data', name='tmpview', sql=self._create_sql(filter), key_column="'IDCELDA'",
                                             workspace='ngd') # las tmp tmb en ngs, supongo
            r.get_response()

    @bcs_session()
    def post(self):
        '''
        {"name":'',
        "wks": '',
        "data":'',
        "file":'',
        "layer_name":'',
        "SRID":'?'} -> layer_name siempre será tmp?
        '''
        db = g.bcs_session.db_session
        t = request.json
        GeographicLayer_Schema= getattr(GeographicLayer,"Schema")()
        geographiclayer = GeographicLayer_Schema.load(t, instance=GeographicLayer())
        db.add(geographiclayer)
        layer_name = geographiclayer.id
        if t.get("data"):
            df = gpd.read_file(t["data"]) # assuming t is a JSON
            self._store_and_publish_gdf(layer_name, df)
        if t.get("file"):
            path = t["file"]
            file_extension = pathlib.Path('my_file.txt').suffix
            if file_extension == ".shp":
                df = gpd.read_file(path)
                self._store_and_publish_gdf(df, layer_name)
            if file_extension in (".tif"):
                r = geoserver_session.create_coveragestore(layer_name=t["id"],
                                         path=path,
                                         workspace=t["wks"])
                r.get_response()
        if t.get("layer_name"):

            layer = geoserver_session.get_layer(layer_name=t["layer_name"])
            # if layer == raster:
            #     r = geoserver_session.create_coveragestore(layer_name=layer_name,
            #                                                path=path,
            #                                                workspace=t["wks"])
            # TODO CAMBIAR FORMATO Y TAL
            self._store_and_publish_gdf(layer_name,layer)

    @bcs_session()
    def put(self, id = None):
        db = g.bcs_session.db_session
        t = request.json
        GeographicLayer_schema = getattr(GeographicRegion, "Schema")()
        issues, geographiclayer, count, status = get_content(db, GeographicLayer, id)
        geographiclayer = GeographicLayer_schema.load(t, instance = geographiclayer)
        geographiclayer.add()
        return ResponseObject(content=geographiclayer,issues = issues, status=status)

    @bcs_session()
    def delete(self, id):
        '''
        delete layer from PostGIS
        ¿and geoserver in case of raster layer?
        @param id:
        @return:
        '''
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        issues, geographiclayer , count, status = get_content(db, GeographicLayer, id)
        if issues == 200:
            layer_name = geographiclayer.id
            # delete layer from geoserver
            r = geoserver_session.delete_layer(layer_name= layer_name, workspace=geographiclayer.wks)
            r.get_response()
            # delete table from postgis
            psqlCursor = pg.cursor()
            dropTableStmt = f"DROP TABLE {layer_name};"
            psqlCursor.execute(dropTableStmt);
            psqlCursor.close();
            # TODO response?
        pg.delete(geographiclayer)
        return ResponseObject() # TODO response


    def _create_sql(self,filter):
        '''
        TODO create sintaxis
        @param filter: filter build by the user in GUI
        @return: sql query for postgis
        '''
        return filter

    def _store_and_publish_gdf(self,layer_name,gdf):
        try:
            gdf.to_postgis(layer_name, postgis_engine,
                          if_exists="fail")  # names have to be uniques or will be stores by id
        except ValueError:
            issues, status = [Issue(IType.INFO, f"layer named as {layer_name} already exists")], 400
            return ResponseObject(issues=issues, status=status, content=None)
        else:
            issues, status = [Issue(IType.INFO, f"layer stored in geo database")], 200
            r = geoserver_session.publish_featurestore(workspace="ngd", store_name='geo_data', pg_table=layer_name)
            r.get_response # todo incluir respuesta de geoserver
            return ResponseObject(issues=issues, status=status, content=None)

register_api(bp_geo, layerAPI, "geo/layers", f"{bcs_api_base}/geo/layers/", pk="layer_id")

class stylesAPI(MethodView):

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