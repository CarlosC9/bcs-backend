import pandas as pd
from flask import Blueprint, request, g
from flask.views import MethodView
from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject, Issue, IType, register_api
from biobarcoding.db_models.geographics import GeographicRegion, Regions
import json
import regex as re

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
        issue = response
        return issue, status
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
        _, status = issues.append(Issue(IType.ERROR, f'no data available')), 200
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
        return ResponseObject(issues= issues,status = status, content= content, content_type= "application/json").get_response()


    @bcs_session()
    def post(self):
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
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
        db.flush()
        return ResponseObject(content=geographicregion, content_type="application/json").get_response()


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
