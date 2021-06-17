import marshmallow.exceptions
import pandas as pd
import geopandas as gpd
from geo.Geoserver import Geoserver
from flask import Blueprint, request, g
from flask.views import MethodView
from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject, Issue, IType, register_api
from biobarcoding.db_models.geographics import GeographicRegion, GeographicLayer, Regions
import json
import pathlib
import regex as re
import requests

from sqlalchemy import create_engine
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT

# postgis_engine = create_engine(app.config['POSTGIS_CONNECTION_STRING'] + "ngd_geoserver")

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
def geoserver_response(response):
    if response:
        status = int(re.findall('[0-9]+', response)[0])
        return Issue(IType.INFO, response.replace(str(status),'')), status
    else:
        return Issue(IType.INFO, "layer succesfully published"), 200

def get_content(session, Feature,issues, id=None):
    content = None
    count = 0
    try:
        content = session.query(Feature)
        if id:
            content = content.filter(Feature.id == id).first()
        else:
            content = content.order_by(Feature.id).all()
            count = len(content)
    except Exception as e:
        print(e)
        _ , status = issues.append(Issue(IType.ERROR, f'READ "{Feature.__name__}" data: The "{Feature.__name__}" data could not be read.')), 500
    if not content:
        _, status = issues.append(Issue(IType.INFO, f'no data available')), 200
    else:
        _, status = issues.append(
            Issue(IType.INFO, f'READ "{Feature.__name__}": The "{Feature.__name__}" data were successfully read')), 200
    return issues,content, count, status

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
        r = ResponseObject()
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        issues = []
        issues, geographicregion, count, status = get_content(db, GeographicRegion,issues, region_id)
        if status == 200 and geographicregion:
            if region_id is None:
                lines = False
                issues, regions, count, status = get_content(pg, Regions,issues, region_id)
            else:
                lines = True
                issues, regions, count, status = get_content(pg, Regions,issues, geographicregion.geo_id)

            bcs_df = pd.read_json(response_to_dataframe(geographicregion),lines = lines)
            postgis_df = pd.read_json(response_to_dataframe(regions), lines = lines)
            bcs_df = bcs_df.set_index(bcs_df["geo_id"]).drop(columns=["geo_id","uuid"])
            postgis_df = postgis_df.set_index(postgis_df["id"]).drop(columns=["uuid","id"])
            content = pd.concat([bcs_df,postgis_df], axis = 1, join="inner")
        else:
            content  = None
        return ResponseObject(issues= issues,status = status, content= content, content_type= "application/json", count = count).get_response()


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
        pg.flush()
        geographicregion = GeographicRegion_schema.load(geographicregion_data, instance = GeographicRegion())
        geographicregion.uuid = regions.uuid
        geographicregion.geo_id = regions.id
        geographicregion.identity_id = g.bcs_session.identity.id
        db.add(geographicregion)
        db.flush()
        return ResponseObject(content= geographicregion, content_type= "application/json").get_response()


    @bcs_session()
    def delete(self,region_id = None):
        issues = []
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        issues, geographicregion, count, status = get_content(db, GeographicRegion, issues, region_id)
        if status == 200 and geographicregion:
            issues, region, count, status = get_content(pg, Regions,issues, geographicregion.geo_id)
            db.delete(geographicregion)
            pg.delete(region)
        return ResponseObject(issues= issues,status = status).get_response()


    @bcs_session()
    def put(self, region_id = None):
        issues = []
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        t = request.json
        GeographicRegion_schema = getattr(GeographicRegion,"Schema")()
        Regions_schema = getattr(Regions, "Schema")()
        issues, geographicregion, count, status = get_content(db, GeographicRegion, issues,  region_id)
        if status == 200 and geographicregion:
            issues, region, count, status = get_content(pg, Regions, issues, geographicregion.geo_id)
            geographicregion = GeographicRegion_schema.load(get_json_from_schema(GeographicRegion, t), instance = geographicregion)
            regions = Regions_schema.load(get_json_from_schema(Regions, t), instance = region)
            db.add(geographicregion)
            pg.add(regions)
        return ResponseObject(issues= issues,status = status, count = count).get_response()


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

    PUT: changing name, attributes  or geometry

    DELETE: delete completely from GeographicLayer table (bcs) postgis and geoserver

    '''
    @bcs_session()
    def get(self, _id = None):
        '''
        - get the list of layers in postgis
        - create a temporal layer from a sql filter
        @param filter: SQL filter in case of temporal layer
            sql = 'SELECT geometry as geom, "DENOMTAX" FROM plantas WHERE "RAREZALOCA" > 1'
            geo.publish_featurestore_sqlview(store_name='geo_data', name='tmp_view3', sql=sql, key_column="DENOMTAX",  (OJO CON LAS COMILLAS)
                                     workspace='ngd')
        @return:
        '''
        # TODO filtrar tmb en el get todo lo que no esté publicado en geoserver?
        from biobarcoding.geo import geoserver_session
        issues = []
        layers = None
        _filter = request.args.get("filter")
        key_col = request.args.get("key_col")
        db = g.bcs_session.db_session
        if _id:
            issues, layers, count, status = get_content(db, GeographicLayer, issues, _id)
            if layers.is_deleted == True:
                layers = None
                _, status = issues.append(Issue(IType.INFO, f'no data available')), 200
        elif not _filter and not key_col:
            issues, layers, count, status = get_content(db, GeographicLayer,issues)
            layers = list(filter(lambda x: (x.is_deleted == False), layers))
        elif _filter and key_col:
            sql = self._create_sql(_filter)
            r = geoserver_session.delete_layer(layer_name='tmpview',workspace='ngd')
            r = geoserver_session.publish_featurestore_sqlview(store_name='geo_data', name='tmpview', sql=sql, key_column=key_col,
                                             workspace='ngd')
            issue, status = geoserver_response(r)
            if not self.check_layer(request):
                _, status = issues.append(Issue(IType.ERROR, f"Error executing request for geoserver")), 500
            issues.append(issue)
            layers = None
        else:
            _, status = issues.append(Issue(IType.ERROR, f"missing data")), 400

        return ResponseObject(issues=issues, status=status, content=layers).get_response()

    @bcs_session()
    def post(self):
        '''
        {"name":"",
        "wks": "",
        "attributes"
        "data":"",
        "file":"",
        "layer_name":'', -> publish tmp view
        " filter ":{"sql":"...","key_col":"..."}
         }, .-> publish view with a name ion geo server and register in bcs as user view
        '''
        issues= []
        db = g.bcs_session.db_session
        t = request.json
        GeographicLayer_Schema= getattr(GeographicLayer,"Schema")()
        geographiclayer_data = get_json_from_schema(GeographicLayer, t)
        geographiclayer = GeographicLayer_Schema.load(geographiclayer_data, instance=GeographicLayer())
        geographiclayer.identity_id = g.bcs_session.identity.id
        db.add(geographiclayer)
        db.flush()
        layer_name = f"layer_{geographiclayer.id}"
        issues, status  = self._post_in_postgis(t, layer_name,issues)
        if status == 200:
            geographiclayer.in_postgis = True
        issues, status, layer_type = self._publish_in_geoserver(t, layer_name, issues)
        if status == 200:
            geographiclayer.published = True
        geographiclayer.layer_type = layer_type
        db.flush()
        return ResponseObject(issues = issues, status=status, content=geographiclayer).get_response()

    @bcs_session()
    def put(self, _id = None):
        issues = []
        db = g.bcs_session.db_session
        t = request.json
        GeographicLayer_schema = getattr(GeographicRegion, "Schema")()
        issues, geographiclayer, count, status = get_content(db, GeographicLayer,issues, _id)
        if geographiclayer:
            try:
                geographiclayer = GeographicLayer_schema.load(t, instance=geographiclayer)
            except marshmallow.exceptions.ValidationError as field:
                geographiclayer_data = get_json_from_schema(GeographicLayer, t)
                geographiclayer = GeographicLayer_schema.load(geographiclayer_data, instance = geographiclayer)
                layer_name = f"layer_{geographiclayer.id}"
                geographiclayer.published = False
                geographiclayer.in_postgis = False
                db.flush()
                if not t.get("wks"):
                    t["wks"] = geographiclayer.wks
                issues, status = self._post_in_postgis(t, layer_name, issues)
                if status == 200:
                    geographiclayer.in_postgis = True
                # re-publish in geoserver is need when the layer is changed (why?)
                # (note that it is not necessary when sql_view layer)
                issues, status, layer_type = self._publish_in_geoserver(t, layer_name, issues)
                if status == 200:
                    geographiclayer.published = True
                    geographiclayer.layer_type = layer_type
                    db.flush()
                else:
                    db.rollback()
        return ResponseObject(content = geographiclayer,issues = issues, status = status).get_response()

    @bcs_session()
    def delete(self, _id):
        '''
        delete layer from PostGIS
        ¿and geoserver in case of raster layer?
        @param id:
        @return:
        '''
        from biobarcoding.geo import geoserver_session
        from biobarcoding import postgis_engine
        issues = []
        db = g.bcs_session.db_session
        issues, geographiclayer , count, status = get_content(db, GeographicLayer,issues, _id)
        if geographiclayer == None:
            return ResponseObject(content=None, issues=issues, status=status).get_response()
        geographiclayer.is_deleted = True
        db.flush()
        if status == 200 and geographiclayer:
            layer_name = f"layer_{str(geographiclayer.id)}"
            # delete layer from geoserver
            r = geoserver_session.delete_layer(layer_name=layer_name, workspace=geographiclayer.wks)
            issue, status = geoserver_response(r)
            issues.append(issue)
            if status ==200:
                geographiclayer.published = False
                db.flush()
            # delete table from postgis
            dropTableStmt = f"DROP TABLE public.{layer_name};"
            try:
                postgis_engine.execute(dropTableStmt)
            except:
                issues.append(Issue(IType.ERROR, f'Error: Error at deleteing {layer_name} layer'))
                status = 404
            else:
                geographiclayer.in_postgis = False
                db.flush()
            if geographiclayer.in_postgis == False and geographiclayer.published == False:
                db.delete(geographiclayer)
        return ResponseObject(content=None, issues = issues, status=status).get_response()


    def _create_sql(self,_filter):
        '''
        TODO create sintaxis
        @param filter: filter build by the user in GUI
        @return: sql query for postgis
        '''
        # key_col = _filter["key"]
        # sql = _filter["sql"]
        return _filter

    def _post_in_postgis(self,t, layer_name, issues):
        from biobarcoding import postgis_engine
        from biobarcoding.geo import geoserver_session
        if t.get("data"):
            df = gpd.GeoDataFrame.from_features(t["data"]['features'])
        elif t.get("layer_name"):
            layer = geoserver_session.get_layer(layer_name=t["layer_name"])
            if isinstance(layer, dict):
                layer = geoserver_session.get_layer(layer_name=t["layer_name"])
                # si no hay tmpview va a lanzar error
                if layer["layer"].get("type") == 'VECTOR':
                    # OPCIÓN 1: pedir los datos directamente desde GEOSERVER mediante WFS (no sé cómo hacerlo)

                    # OPCIÓN 2: Hacer el query a postgis y guardarlo en postgis y geoserver
                    response = requests.get(layer["layer"]["resource"]["href"], auth=(geoserver_session.username, geoserver_session.password))
                    if response.status_code == 200:
                        tmpview_json = response.text
                        tmpview_dict = json.loads(tmpview_json)
                        sql = tmpview_dict["featureType"]["metadata"]["entry"]["virtualTable"]["sql"]
                        key_col = tmpview_dict["featureType"]["metadata"]["entry"]["virtualTable"]["keyColumn"]
                        df = gpd.GeoDataFrame.from_postgis(sql, postgis_engine)
                    else:
                        _, status = issues.append(Issue(IType.INFO, f"no tmpview available")), response.status_code
                        return issues,status
            else:  # go to geoserver
                return issues, None
        elif t.get("file"):
            path = t["file"]
            file_extension = pathlib.Path(path).suffix
            if file_extension == ".shp":
                df = gpd.read_file(path)

            else: # go to geoserver
                return issues, None
        else:
            return issues, None
        try:
            df.to_postgis(layer_name, postgis_engine, if_exists="replace")
        except ValueError as e:
            _, status =  issues.append(Issue(IType.INFO, f"layer named as {layer_name} error {e}")), 400
        else:
            _, status =  issues.append(Issue(IType.INFO, f"layer stored in geo database")), 200
        return issues, status


    def _publish_in_geoserver(self, t, layer_name, issues):
        from biobarcoding.geo import geoserver_session
        layer_type = None
        # TODO REFACTOR
        layer_type = "vector"
        if t.get("data"):
            r = geoserver_session.publish_featurestore(workspace=t['wks'], store_name='geo_data',
                                                       pg_table=layer_name)
        elif t.get("layer_name"):
            layer = geoserver_session.get_layer(layer_name=t["layer_name"])
            if isinstance(layer,dict):
                if layer["layer"].get("type") == 'VECTOR':
                    r = geoserver_session.publish_featurestore(workspace=t['wks'], store_name='geo_data',
                                                               pg_table=layer_name)
            else:
                r = "no tmpview available 500 "

        elif t.get("filter"):
            sql = t["filter"]
            key_col = t["key_col"]
            r = geoserver_session.publish_featurestore_sqlview(store_name='geo_data', name=layer_name, sql=sql,
                                                               key_column=key_col,
                                                               workspace=t['wks'])
            layer_type = "sqlview"

        elif t.get('file'):
            path = t.get('file')
            file_extension = pathlib.Path(path).suffix
            if file_extension in (".tif"):
                r = geoserver_session.create_coveragestore(layer_name=layer_name,
                                                           path=path,
                                                           workspace=t["wks"])
                layer_type = "raster"
            elif file_extension in (".shp"):
                r = geoserver_session.publish_featurestore(workspace=t['wks'],
                                                           store_name='geo_data',
                                                           pg_table=layer_name)
            else:
                r = f'READ "no valid file extension: {file_extension} 400'
        else:
            r = "No  valid geographic data specified, 400"

        issue, status = geoserver_response(r)
        issues.append(issue)
        return issues, status, layer_type

    def check_layer(self, request):
        from biobarcoding.geo import geoserver_session
        layer = geoserver_session.get_layer(layer_name = "tmpview")
        tmpview_json = requests.get(layer["layer"]["resource"]["href"],
                                    auth=(geoserver_session.username, geoserver_session.password)).text
        tmpview_dict = json.loads(tmpview_json)
        sql = tmpview_dict["featureType"]["metadata"]["entry"]["virtualTable"]["sql"].rstrip("\n")
        key_col = tmpview_dict["featureType"]["metadata"]["entry"]["virtualTable"]["keyColumn"]
        if sql == request.args.get("filter") and key_col == request.args.get("key_col"):
            return True
        else:
            return False

    def _publish_view(self, t, layer_name):
        pass

    def _post_view(self, t, layer_name):
        pass





view_func = layerAPI.as_view("geo/layers")
bp_geo.add_url_rule(f"{bcs_api_base}/geo/layers/", defaults={"_id" : None}, view_func=view_func, methods=['GET'])
bp_geo.add_url_rule(f"{bcs_api_base}/geo/layers/<int:_id>", view_func=view_func, methods=['GET'])
bp_geo.add_url_rule(f"{bcs_api_base}/geo/layers/", view_func=view_func, methods=['POST'])
bp_geo.add_url_rule(f'{bcs_api_base}/geo/layers/<int:_id>', view_func=view_func, methods=['PUT', 'DELETE'])

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