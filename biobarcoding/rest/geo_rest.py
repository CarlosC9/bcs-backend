import io
import json
import os
import pathlib
import tempfile
from typing import List, Tuple
from urllib.parse import urlparse

import geopandas as gpd
import pandas as pd
import regex as re
import requests
import seaborn as sns
from flask import Blueprint, request, g, Response
from flask.views import MethodView
from marshmallow import EXCLUDE
from matplotlib.colors import rgb2hex

from ..authentication import n_session
from ..db_models.geographics import GeographicRegion, Regions, GeographicLayer
from ..geo import workspace_names, postgis_store_name
from ..geo.biota import read_biota_file, generate_pda_species_file_from_layer, import_pda_result
from . import app_api_base, ResponseObject, Issue, IType, register_api, app_proxy_base, filter_parse
from ..services.files import get_file_contents

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


def geoserver_response(response) -> Tuple[Issue, int]:
    """
    Encode GeoServer response into an Issue

    :param response:
    :return:
    """
    if response:
        if response == "Could not acquire reader for coverage":
            status = 400
        else:
            status = int(re.findall('[0-9]+', response)[0])
        return Issue(IType.INFO, response.replace(str(status), '')), status
    else:
        return Issue(IType.INFO, "layer succesfully published"), 200


def get_content(session, feature_class, issues, id_=None, filter_=None):
    def __aux_own_filter(filt_):
        """
            Example clause:
                {"tags": "biota"} meaning "the layer contains 'biota' as one of the tags"
        """
        clause = []

        if 'tags' in filt_:
            v = filt_.get('tags')
            # 3 implementations
            # from sqlalchemy import text
            # clause.append(text(f"attributes->'tags' ? :n").params(n=v))
            # clause.append(GeographicLayer.attributes["tags"].op("?")(v))
            clause.append(GeographicLayer.attributes["tags"].has_key(v))

        return clause

    content = None
    count = 0
    try:
        content = session.query(feature_class)
        if id_:
            content = content.filter(feature_class.id == id_).first()
        else:
            if filter_:
                content = content.filter(filter_parse(GeographicLayer, filter_, __aux_own_filter))
            content = content.order_by(feature_class.id).all()
            count = len(content)
    except Exception as e:
        print(e)
        _, status = issues.append(Issue(IType.ERROR,
                                        f'READ "{feature_class.__name__}" data: The "{feature_class.__name__}"'
                                        f' data could not be read.')), 500
    if not content:
        _, status = issues.append(Issue(IType.INFO, f'no data available')), 200
    else:
        _, status = issues.append(
            Issue(IType.INFO, f'READ "{feature_class.__name__}": The "{feature_class.__name__}"'
                              f' data were successfully read')), 200
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


def get_json_from_schema(entity, input_):
    entity_schema = getattr(entity, "Schema")()
    t_json = entity_schema.dumps(input_)
    return json.loads(t_json)


def generate_ramp_sld_file(
        style_name: str,
        column_name: str,
        values: List[float],
        n_intervals: int = 5,
        color_ramp: str = None,
        geom_type: str = "polygon",
):
    """


    Base on function "Style.classified_xml"

    :param style_name:
    :param column_name:
    :param values:
    :param n_intervals:
    :param color_ramp:
    :param geom_type:
    :return:
    """
    max_value = max(values)
    min_value = min(values)
    diff = max_value - min_value
    n = n_intervals
    interval = diff / 5
    palette = sns.color_palette(color_ramp, int(n))
    palette_hex = [rgb2hex(i) for i in palette]
    # interval = N/4
    # color_values = [{value: color} for value, color in zip(values, palette_hex)]
    # print(color_values)
    rule = ""
    for i, color in enumerate(palette_hex):
        rule += """
            <sld:Rule>
                <sld:Name>{1}</sld:Name>
                <sld:Title>{4}</sld:Title>
                <ogc:Filter xmlns:ogc="http://www.opengis.net/ogc">
                    <ogc:And>
                    <ogc:PropertyIsGreaterThan{6}>
                        <ogc:PropertyName>{0}</ogc:PropertyName>
                        <ogc:Literal>{4}</ogc:Literal>
                    </ogc:PropertyIsGreaterThan{6}>
                    <ogc:PropertyIsLessThanOrEqualTo>
                        <ogc:PropertyName>{0}</ogc:PropertyName>
                        <ogc:Literal>{5}</ogc:Literal>
                    </ogc:PropertyIsLessThanOrEqualTo>
                    </ogc:And>
                </ogc:Filter>
                <sld:PolygonSymbolizer>
                    <sld:Fill>
                      <sld:CssParameter name="fill">{2}</sld:CssParameter>
                    </sld:Fill>
                    <sld:Stroke>
                      <sld:CssParameter name="stroke">{3}</sld:CssParameter>
                      <sld:CssParameter name="stroke-width">1</sld:CssParameter>
                      <sld:CssParameter name="stroke-linejoin">mitre</sld:CssParameter>
                    </sld:Stroke>
                </sld:PolygonSymbolizer>
            </sld:Rule>

        """.format(
            column_name,
            style_name,
            color,
            "#000000",
            min_value + interval * i,
            min_value + interval * (i + 1),
            "OrEqualTo" if i == 0 else ""
        )

    style = """
            <sld:StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:sld="http://www.opengis.net/sld" xmlns:gml="http://www.opengis.net/gml" xmlns:ogc="http://www.opengis.net/ogc" version="1.0.0">
                <sld:NamedLayer>
                    <sld:Name>{0}</sld:Name>
                    <sld:UserStyle>
                    <sld:Name>{0}</sld:Name>
                        <sld:FeatureTypeStyle>
                            {1}
                        </sld:FeatureTypeStyle>
                    </sld:UserStyle>
                </sld:NamedLayer>
            </sld:StyledLayerDescriptor>
        """.format(
        style_name, rule
    )

    with open("style.sld", "w") as f:
        f.write(style)


def create_and_publish_ramp_style(gs_session,
                                  wkspc,
                                  layer_name, attribute, min_value, max_value,
                                  style_name, color_ramp: str = "RdYlGn_r", cmap_type: str = "ramp",
                                  number_of_classes: int = 5, overwrite: bool = False):
    sld_version = "1.0.0"
    generate_ramp_sld_file(style_name=style_name,
                           column_name=attribute,
                           values=[min_value, max_value],
                           n_intervals=number_of_classes,
                           color_ramp=color_ramp)
    # Style.classified_xml(style_name=style_name,
    #                      column_name=attribute,
    #                      values=[min_value, max_value],
    #                      color_ramp=color_ramp)
    # sld_version = "1.1.0"

    gs_session.upload_style("style.sld", style_name, wkspc, sld_version=sld_version, overwrite=overwrite)
    os.remove("style.sld")
    gs_session.publish_style(layer_name, style_name, wkspc)


class LayersAPI(MethodView):
    """
    GET:    the list of layers info from Geographiclayer table on bcs
            publish a result of filtering a layer (with geoserver_session.publish_featurestore_sqlview(sql_string))
            TODO: Move the dynamic layer behavior to PUT, because it is a -assumed preexisting- virtual layer
                  that is modified. PUT /api/geo/layers/0?definition=... ("0" would be a special ID)

    POST:   Three post cases:
            1. post all new vector layer from geoJSON format (PDA case)
            2. Import a Raster o Vector layer from a file (via command)
            3. Save a temporal layer created by geoserver_session.publish_featurestore_sqlview(sql_string)
            create and store in GeographicLayer table as item
            then in postgis as table (format)
            and then publish  in geoserver as new layer

    PUT:    changing name, attributes or geometry

    DELETE: delete completely from GeographicLayer table (bcs) postgis and geoserver
    """

    # Some defaults
    kwargs = {}
    issues = []
    status = int()

    def _export(self, sess, lay: GeographicLayer, _format):
        if _format not in ["nexus", "pda_simple", "species", "species_canon"]:
            self.issues.append(
                Issue(IType.ERROR, f'Could not export Biota layer {lay.name} (internal name "{lay.geoserver_name}") '
                                   f'to format "{_format}". Supported "pda_simple", "nexus", "species", "species_canon"',
                      f"Export {lay.name}"))
            return None
        # TODO Check "tags" attribute, the layer should have a "Biota" tag
        _ = generate_pda_species_file_from_layer(sess, lay.id, lay.geoserver_name, _format)
        if _ is None:
            self.issues.append(
                Issue(IType.ERROR, f'Could not export Biota layer {lay.name} (internal name "{lay.geoserver_name}").',
                      f"Export {lay.name} as Simple PDA text file"))
            return None
        elif _ == "":
            self.issues.append(Issue(IType.ERROR,
                                     f'Empty Biota layer {lay.name} (internal name "{lay.geoserver_name}").',
                                     f"Export {lay.name}"))
            return None
        else:
            return _

    @n_session()
    def get(self, _id=None, _format=None):
        """
        - get the list of layers in postgis
        - create a temporal layer from a sql filter
        if query params "filter" and "key_column": SQL filter, to calculate a temporary layer
            geo.publish_featurestore_sqlview(
                store_name='geo_data',
                name='tmp_view3',
                sql='SELECT geometry as geom, "DENOMTAX" FROM plantas WHERE "RAREZALOCA" > 1',
                key_column="DENOMTAX",  (OJO CON LAS COMILLAS)
                workspace='ngd')

        USAGE (CURL EXAMPLE):
        export API_BASE_URL=http://localhost:5000/api
        curl --cookie-jar bcs-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        curl --cookie bcs-cookies.txt "$API_BASE_URL/geo/layers/"
        curl --cookie bcs-cookies.txt "$API_BASE_URL/geo/layers/?filter=%7B%22tags%22%3A%20%22biota%22%7D"

        curl --cookie bcs-cookies.txt "$API_BASE_URL/geo/layers/1"

        Specially, a Biota layer can be exported as one of 3 supported formats. The first two are for PDA, the third
        just enumerates available species:

        curl --cookie bcs-cookies.txt "$API_BASE_URL/geo/layers/1.pda_simple"
        curl --cookie bcs-cookies.txt "$API_BASE_URL/geo/layers/1.nexus"
        curl --cookie bcs-cookies.txt "$API_BASE_URL/geo/layers/1.species"
        curl --cookie bcs-cookies.txt "$API_BASE_URL/geo/layers/1.species_canon"

        @param _id: ID of a layer; empty for ALL layers (if no query params specified)
        @param _format: export format
        @return:
        """
        from biobarcoding.geo import geoserver_session
        self.issues = []
        layer = None
        _filter = request.args.get("filter")
        if _filter:
            _filter = json.loads(_filter)
        key_col = request.args.get("key_col")
        db_sess = g.n_session.db_session
        if _id:  # A layer
            self.issues, layer, count, self.status = get_content(db_sess, GeographicLayer, self.issues, _id)
            if layer and layer.is_deleted:
                layer = None
                _, self.status = self.issues.append(Issue(IType.INFO, f'no data available')), 200
            else:
                if _format:
                    content = self._export(db_sess, layer, _format=_format)
                    if content:
                        return Response(content, mimetype=f"text/{_format}", status=200)
                else:
                    if layer:
                        serializer = layer.Schema()
                        serializer.dump(layer)
        elif not key_col:  # ALL LAYERS layers (maybe filtered)
            self.issues, layer, count, self.status = get_content(db_sess, GeographicLayer, self.issues, filter_=_filter)
            if layer:
                layer = list(filter(lambda x: (x.is_deleted is False), layer))
        elif _filter and key_col:  # Temporary layer
            # TODO Ensure Temporary layer is also registered as BCS GeographicLayer
            tmp_view = f'tmpview_{g.n_session.identity.id}'
            sql = self._create_sql(_filter)
            r = geoserver_session.delete_layer(layer_name=tmp_view, workspace=workspace_names[1])
            print(r)
            r = geoserver_session.publish_featurestore_sqlview(store_name=postgis_store_name, name=tmp_view, sql=sql,
                                                               key_column=key_col,
                                                               workspace=workspace_names[1])
            print(r)
            issue, self.status = geoserver_response(r)
            if not self._check_layer():
                _, self.status = self.issues.append(Issue(IType.ERROR, f"Error executing request for geoserver")), 500
            self.issues.append(issue)
            layer = dict(tmpview=tmp_view, sql=sql, key_col=key_col)
        else:  # No information to elaborate a response
            _, self.status = self.issues.append(Issue(IType.ERROR, f"missing data")), 400

        return ResponseObject(issues=self.issues, status=self.status, content=layer).get_response()

    @staticmethod
    def _exclude(c):
        """
        Check if column "c" has to be excluded from the list of attributes of a layer showed to users

        :param c:
        :return:
        """
        s = c.lower()
        return s.startswith("id") or s.endswith("id") or "codigo" in s or s in ["coordx", "coordy", "geom", "geometry"]

    def _create_properties_and_geoserver_styles(self, gdf, wks, layer_name, lc_attributes):
        from biobarcoding.geo import geoserver_session
        _ = []
        for c in gdf.columns:
            if self._exclude(c):
                continue
            if lc_attributes:
                c = c.lower()
            tmp = gdf[c].values
            try:
                min_v = float(min(tmp))
                max_v = float(max(tmp))
                p_type = "numeric"
                style_name = f"{layer_name}_{c}"
            except:
                min_v = None
                max_v = None
                p_type = "string"
                style_name = None

            # Create style for properties (_publish_in_geoserver
            if style_name:
                create_and_publish_ramp_style(geoserver_session,
                                              wkspc=wks,
                                              layer_name=layer_name,
                                              attribute=c,
                                              min_value=min_v, max_value=max_v,
                                              number_of_classes=7,
                                              style_name=style_name)

            _.append(dict(name=c, type=p_type, style=style_name, min=min_v, max=max_v))
        return _

    @n_session()
    def post(self):
        """
        1. import a raster file (tif and twf file)
        2. import vector layer from SHP file (import shp shx...)
        3. import layer from GeoJSON data
        4. create a new layer from a temp_view
        5. create a new vector file from shape file and convert data to json
        {"name":"",
        "wks": "ngd",
        "attributes" -> json attributes like {"tags"["tag1", "tag2"]
        "data":"", -> import layer from geojson data
        "layer_name":'', -> publish tmp view
        "convert_to: "geojson", -> convert Biota file to geojson (very long process)
         }

        USAGE (CURL EXAMPLE):
          1) Define base URL, 2) Login, 3) Post multipart-form with a file and a JSON to "/api/geo/layers"

        export API_BASE_URL=http://localhost:5000/api
        curl --cookie-jar bcs-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        curl --cookie bcs-cookies.txt -F "layer_file=@/home/rnebot/GoogleDrive/AA_NEXTGENDEM/plantae_canarias/Plantas.zip;type=application/zip" -F "metadata={\"name\": \"capa_1\", \"wks\": \"ngd\", \"attributes\": {\"tags\": [\"tag1\", \"tag2\"]}};type=application/json" "$API_BASE_URL/geo/layers/"

        POST from FilesAPI:
        First, upload a file to FilesAPI:
        curl --cookie bcs-cookies.txt -H "Content-Type: application/zip" -XPUT --data-binary @/home/rnebot/GoogleDrive/AA_NEXTGENDEM/plantae_canarias/Plantas.zip "$API_BASE_URL/files/f1/f2/f3/plantas.zip.content"
        Then, do the POST, like
        curl --cookie bcs-cookies.txt -XPOST "$API_BASE_URL/geo/layers/?filesAPI=%2Ff1%2Ff2%2Ff3%2Fplantas.zip&job_id=6"
        """
        db = g.n_session.db_session
        self.issues = []
        self.kwargs = self._build_args_from_request_data()
        system_layer = True  # or user layer
        self.kwargs["wks"] = workspace_names[0] if system_layer else workspace_names[1]
        # Put self.kwargs into a geographic layer
        geographic_layer_schema = getattr(GeographicLayer, "Schema")()
        geographic_layer_data = get_json_from_schema(GeographicLayer, self.kwargs)
        # Create object in memory
        geographic_layer = geographic_layer_schema.load(geographic_layer_data, instance=GeographicLayer())
        geographic_layer.identity_id = g.n_session.identity.id
        # Persist object in BCS DB
        db.add(geographic_layer)
        db.flush()
        # Receive file
        if self._file_posted():
            # Set URL
            tmp = urlparse(request.base_url)
            base_url = f"{tmp.scheme}://{tmp.netloc}"
            geographic_layer.wms_url = f"{base_url}{app_proxy_base}/geoserver/wms"
            # Store file locally
            self._receive_and_prepare_file()
            # Store file in PostGIS
            lower_case_attributes = True
            layer_name = f"layer_{geographic_layer.id}"  # (internal) Geoserver layer name
            status, gdf = self._post_in_postgis(layer_name, lower_case_attributes)
            # Delete temporary file
            os.remove(self.kwargs["path"])
            # Publish layer in Geoserver
            if status == 200:
                geographic_layer.in_postgis = True

            status, layer_type = self._publish_in_geoserver(layer_name)
            # if self.kwargs.get("property"):
            #     prop = self.kwargs["property"]
            #     if lower_case_attributes:
            #         prop = prop.lower()
            #     tmp = gdf[prop].values
            #     style_name = f"{layer_name}_{prop}"
            #     status, layer_type = self._publish_in_geoserver(layer_name, prop, style_name, min(tmp), max(tmp))
            # else:
            #     prop = None
            #     status, layer_type = self._publish_in_geoserver(layer_name)

            if status == 200:
                geographic_layer.published = True
                geographic_layer.geoserver_name = layer_name
                geographic_layer.properties = self._create_properties_and_geoserver_styles(gdf,
                                                                                           self.kwargs["wks"],
                                                                                           layer_name,
                                                                                           lower_case_attributes)
            geographic_layer.layer_type = layer_type
            db.flush()
        return ResponseObject(issues=self.issues, status=self.status, content=geographic_layer).get_response()

    @n_session()
    def put(self, _id=None):
        """
        Update geographic layer metadata

        {"name":"",
        "wks": "ngd",
        "attributes" -> json attributes like {"tags"["tag1", "tag2"]
        "data":"", -> import layer from geojson data
        "layer_name":'', -> publish tmp view
         }

        USAGE (CURL EXAMPLE):
          1) Define base URL, 2) Login, 3) Post multipart-form with a file and a JSON to "/api/geo/layers"

        export API_BASE_URL=http://localhost:5000/api
        curl --cookie-jar bcs-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        curl --cookie bcs-cookies.txt -X PUT -F "metadata={\"name\": \"capa_12\", \"wks\": \"ngd\", \"attributes\": {\"tags\": [\"tag1\", \"tag2\", \"tag3\"]}};type=application/json" "$API_BASE_URL/geo/layers/"

        :param _id:
        :return:
        """
        db = g.n_session.db_session
        self.issues = []
        self.kwargs = self._build_args_from_request_data()
        system_layer = True
        self.kwargs["wks"] = workspace_names[0] if system_layer else workspace_names[1]
        geographic_layer_schema = getattr(GeographicLayer, "Schema")()
        geographic_layer_data = get_json_from_schema(GeographicLayer, self.kwargs)
        # Remove unmodifiable keys
        for k in ["geoserver_name", "id", "identity_id", "in_postgis", "is_deleted", "layer_type", "published", "uuid",
                  "wks"]:
            if k in geographic_layer_data:
                del geographic_layer_data[k]
        self.issues, geographic_layer, count, self.status = get_content(db, GeographicLayer, self.issues, _id)
        if geographic_layer:
            # Update layer metadata
            geographic_layer = geographic_layer_schema.load(geographic_layer_data,
                                                            instance=geographic_layer,
                                                            partial=True, unknown=EXCLUDE)
            # TODO Update layer content disabled (not defined how to
            # if self.kwargs.get("data") or request.files:
            #     layer_name = f"layer_{geographic_layer.id}"
            #     geographic_layer.published = False
            #     geographic_layer.in_postgis = False
            #     db.flush()
            #     if request.files:
            #         self._receive_and_prepare_file()
            #     if not self.kwargs.get("wks"):
            #         self.kwargs["wks"] = geographic_layer.wks
            #     status, gdf = self._post_in_postgis(layer_name)
            #     if status == 200:
            #         geographic_layer.in_postgis = True
            #     # re-publish in geoserver is need when the layer is changed (why?)
            #     # (note that it is not necessary when sql_view layer)
            #     if self.kwargs.get("property"):
            #         prop = self.kwargs["property"]
            #         tmp = gdf[prop].values
            #         status, layer_type = self._publish_in_geoserver(layer_name, prop, min(tmp), max(tmp))
            #     else:
            #         status, layer_type = self._publish_in_geoserver(layer_name)
            #     if status == 200:
            #         geographic_layer.published = True
            #         geographic_layer.layer_type = layer_type
            #         db.flush()
            #     else:
            #         db.rollback()
        else:
            if _id is None:
                self.issues.append(Issue(IType.ERROR, '<id> not specified'))
            else:
                self.issues.append(Issue(IType.ERROR, f'Layer with specified id ({_id}) was not found'))
        return ResponseObject(content=geographic_layer, issues=self.issues, status=self.status).get_response()

    @n_session()
    def delete(self, _id):
        """
        delete layer from PostGIS
        ¿and geoserver in case of raster layer?
        @param _id:
        @return:
        """
        from biobarcoding import postgis_engine
        from biobarcoding.geo import geoserver_session
        self.issues = []
        db = g.n_session.db_session
        self.issues, geographic_layer, count, self.status = get_content(db, GeographicLayer, self.issues, _id)
        if geographic_layer is None:
            return ResponseObject(content=None, issues=self.issues, status=self.status).get_response()
        geographic_layer.is_deleted = True
        db.flush()
        if self.issues == 200 and geographic_layer:
            layer_name = f"layer_{str(geographic_layer.id)}"
            # delete layer from geoserver
            r = geoserver_session.delete_layer(layer_name=layer_name, workspace=geographic_layer.wks)
            issue, self.status = geoserver_response(r)
            self.issues.append(issue)
            if self.issues == 200:
                geographic_layer.published = False
                db.flush()
            # delete table from postgis
            drop_table_stmt = f"DROP TABLE public.{layer_name};"
            try:
                postgis_engine.execute(drop_table_stmt)
            except:
                self.issues.append(Issue(IType.ERROR, f'Error: Error at deleting {layer_name} layer'))
                self.status = 404
            else:
                geographic_layer.in_postgis = False
                db.flush()
            if geographic_layer.in_postgis == False and geographic_layer.published == False:
                db.delete(geographic_layer)
        return ResponseObject(content=None, issues=self.issues, status=self.status).get_response()

    def _create_sql(self, _filter):
        """
        TODO create sintaxis
        @param _filter: filter build by the user in GUI
        @return: sql query for postgis
        """
        # key_col = _filter["key"]
        # sql = _filter["sql"]
        # _filter = _filter.encode('utf-8')
        return _filter

    def _post_in_postgis(self, layer_name, lower_columns=True):
        """
        Create (store) feature layer in PostGIS
        NOTE: implicit parameters in self.kwargs

        :param layer_name:
        :param lower_columns: rename all columns to lower case
        :return: None if error, status if success
        """
        from biobarcoding import postgis_engine
        from biobarcoding.geo import geoserver_session

        if self.kwargs["path"] != "":
            df = self._read_vector_file()
            if not isinstance(df, gpd.GeoDataFrame):
                return None
        elif "data" in self.kwargs.keys():
            df = gpd.GeoDataFrame.from_features(self.kwargs["data"]['features'])
        elif "layer_name" in self.kwargs.keys():
            layer = geoserver_session.get_layer(layer_name=self.kwargs["layer_name"])
            if isinstance(layer, dict):
                layer = geoserver_session.get_layer(layer_name=self.kwargs["layer_name"])
                # si no hay tmpview va a lanzar error
                if layer["layer"].get("type") == 'VECTOR':
                    # TODO OPTION 1: export data and execute a Python function
                    # OPTION 2: execute a SQL view statement into PostGIS
                    df = self._get_df_from_view()
            else:  # go to geoserver
                return None
        try:
            # Make PK for sqlview work properly using id_column as PK
            idx_column = self.kwargs.get("idx_field")
            if idx_column:
                if idx_column in df.columns:
                    # df = df.drop(columns=idx_column, axis=1)
                    del df[idx_column]
                df.to_postgis(layer_name, postgis_engine, if_exists="replace", index=True, index_label=idx_column)
                with postgis_engine.connect() as con:
                    con.execute(f'ALTER TABLE {layer_name} ADD PRIMARY KEY ({idx_column});')
            else:
                # All columns to lower case
                df.columns = [c.lower() for c in df.columns]
                # Write to PostGIS
                df.to_postgis(layer_name, postgis_engine, if_exists="replace")
        except ValueError as e:
            _, self.status = self.issues.append(Issue(IType.INFO, f"Layer named as {layer_name} error {e}")), 400
        else:
            self.kwargs[
                "data"] = "dummy"  # "add geoserver layer" function looks for any value in "data" to create a layer that is stored in PostGIS
            _, self.status = self.issues.append(Issue(IType.INFO, f"layer stored in geo database")), 200
        return self.status, df

    def _publish_in_geoserver(self, layer_name, attribute=None, style_name=None, min_value=None, max_value=None):
        """
        Publish a layer -already in a registered store- in Geoserver
        NOTE: implicit parameters in self.kwargs

        NOTE!!!: min_value and max_value were previously a "range" list with two elements. For some reason,
                 using [min(tmp), max(tmp)] to prepare the parameter in the call, breaks "list" reserved word

        :param layer_name:
        :param attribute: attribute name for default style
        :param style_name: name of the style to create
        :param min_value: mininum for default style
        :param max_value: maximum for default style
        :return:
        """
        from biobarcoding.geo import geoserver_session
        layer_type = "vector"
        if "data" in self.kwargs.keys():
            r = geoserver_session.delete_layer(layer_name=layer_name, workspace=self.kwargs['wks'])
            print(r)
            r = geoserver_session.publish_featurestore(workspace=self.kwargs['wks'],
                                                       store_name=postgis_store_name,
                                                       pg_table=layer_name)
            if attribute:
                create_and_publish_ramp_style(geoserver_session,
                                              wkspc=self.kwargs['wks'],
                                              layer_name=layer_name,
                                              attribute=attribute,
                                              min_value=min_value, max_value=max_value,
                                              number_of_classes=7,
                                              style_name=style_name)
        elif "layer_name" in self.kwargs.keys():
            layer = geoserver_session.get_layer(layer_name=self.kwargs["layer_name"])
            if isinstance(layer, dict):
                if layer["layer"].get("type") == 'VECTOR':
                    r = geoserver_session.publish_featurestore(workspace=self.kwargs['wks'],
                                                               store_name=postgis_store_name,
                                                               pg_table=layer_name)
            else:
                r = "no view available 500 "

        elif self.kwargs.get("filter"):
            sql = self.kwargs["filter"]
            key_col = self.kwargs["key_col"]
            r = geoserver_session.publish_featurestore_sqlview(store_name=postgis_store_name,
                                                               name=layer_name,
                                                               sql=sql,
                                                               key_column=key_col,
                                                               workspace=self.kwargs['wks'])
            layer_type = "sqlview"
        elif request.files:
            r, layer_type = self._read_raster_file(layer_name)
        else:
            r = "No  valid geographic data specified, 400"

        issue, self.status = geoserver_response(r)
        self.issues.append(issue)
        return self.status, layer_type

    def _tmpview_sql(self, layer_name):
        from biobarcoding.geo import geoserver_session
        layer = geoserver_session.get_layer(layer_name=layer_name)
        if not isinstance(layer, dict):
            return None, None
        response = requests.get(layer["layer"]["resource"]["href"],
                                auth=(geoserver_session.username, geoserver_session.password))
        if response.status_code == 200:
            tmp_view_json = response.text
            tmp_view_dict = json.loads(tmp_view_json)
            sql = tmp_view_dict["featureType"]["metadata"]["entry"]["virtualTable"]["sql"].rstrip("\n")
            key_col = tmp_view_dict["featureType"]["metadata"]["entry"]["virtualTable"]["keyColumn"]
            return sql, key_col
        else:
            return None, None

    def _check_layer(self):
        sql, key_col = self._tmpview_sql("tmpview")
        if sql == request.args.get("filter") and key_col == request.args.get("key_col"):
            return True
        else:
            return False

    def _get_df_from_view(self):
        sql, key_col = self._tmpview_sql(self.kwargs["layer_name"])
        if sql:
            from biobarcoding import postgis_engine
            df = gpd.GeoDataFrame.from_postgis(sql, postgis_engine)
            return df
        else:
            _, self.status = self.issues.append(Issue(IType.INFO, f"no tmpview available")), 400
            return None

    @staticmethod
    def _build_args_from_request_data():
        args = {}
        if request.json:
            args.update(request.json)
        elif request.form:
            if request.form.get("metadata"):
                args.update(json.loads(request.form["metadata"]))
        elif request.values:
            t = request.values.to_dict()
            for key, item in t.items():
                if isinstance(item, str):
                    try:
                        item = json.loads(item)
                        args.update(item)
                    except json.decoder.JSONDecodeError:
                        args.update(t)
                        break
        if "attributes" not in args:
            # Define default attributes
            args["attributes"] = dict(tags=[])
        else:
            if "tags" not in args["attributes"]:
                args["attributes"]["tags"] = []
        if "filesAPI" in args:
            args["attributes"]["job_file"] = args["filesAPI"]
        if "job_id" in args:
            args["attributes"]["job_id"] = args["job_id"]
            del args["job_id"]
        if "name" not in args:
            # Define default name
            if "job_id" in args["attributes"]:
                args["name"] = f'Capa salida del job {args["attributes"]["job_id"]}'
        return args

    def _file_posted(self):
        return len(request.files) > 0 or "filesAPI" in self.kwargs

    def _receive_and_prepare_file(self):
        import os
        formats = [".tif", ".gpkg", ".zip", ".json", ".geojson"]
        path = ""
        folder = tempfile.mkdtemp(prefix="bcs_")
        from werkzeug.utils import secure_filename

        # Download to a local file if self.kwargs.data contains a FilesAPI path
        if "filesAPI" in self.kwargs:
            c_type, contents = get_file_contents(g.n_session.db_session, self.kwargs["filesAPI"])
            from mimetypes import guess_extension
            temp_name = tempfile.NamedTemporaryFile(dir=folder, delete=False)
            file_path = os.path.join(folder, f"{temp_name.name}{guess_extension(c_type)}")
            with open(file_path, "wb") as f:
                f.write(contents)
            path = file_path
        elif len(request.files) > 0:
            files = request.files.to_dict(flat=False)
            for _, value in files.items():  # TODO It will consider just the last file
                for file in value:
                    file_path = os.path.join(folder, secure_filename(file.filename))
                    file.save(file_path)
                    if pathlib.Path(file_path).suffix in formats:
                        path = file_path
        elif len(request.data) > 0:
            if request.content_type == "application/zip":  # Shapefile
                file = f"submitted_file.zip"
            elif request.content_type == "application/geopackage":
                file = f"submitted_file.gpkg"
            file_path = os.path.join(folder, secure_filename(file))
            buffer = bytes(io.BytesIO(request.get_data()).getbuffer())
            with open(file_path, "wb") as f:
                f.write(buffer)
            path = file_path

        if path == "":
            self.status = 400
            self.issues.append(Issue(IType.ERROR, 'unrecognised input file'))
        self.kwargs["path"] = path

    def _read_vector_file(self):
        path = self.kwargs["path"]
        if self.kwargs.get("convert_to"):  # TODO RNEBOT: remove, too specific. Instead, process layer BEFORE submission
            if self.kwargs["convert_to"] == "geojson":
                gdf = read_biota_file(path)
                return gdf
            else:
                self.issues.append(Issue(IType.ERROR, "No valid data to convert"))
                self.status = 400
                return None
        else:  # Main
            try:
                gdf = gpd.read_file(path)
                return gdf
            except:
                try:
                    gdf = import_pda_result(path, session=g.n_session.db_session)
                    return gdf
                except:
                    return None

    def _read_raster_file(self, layer_name):
        from biobarcoding.geo import geoserver_session
        path = self.kwargs["path"]
        file_extension = pathlib.Path(path).suffix
        if file_extension in (".tif"):
            self.kwargs["wks"] = 'ngd'
            r = geoserver_session.create_coveragestore(layer_name=layer_name,
                                                       path=path,
                                                       workspace=self.kwargs["wks"])
            print(r)
            r = geoserver_session.create_coveragestyle(raster_path=path,
                                                       style_name=layer_name + '_style',
                                                       workspace=self.kwargs["wks"])
            print(r)
            r = geoserver_session.publish_style(layer_name=layer_name, style_name=layer_name + '_style',
                                                workspace=self.kwargs["wks"])
            print(r)
            # TODO improve error control when importing to geoserver
            layer_type = "raster"
        elif file_extension in (".shp"):
            r = geoserver_session.publish_featurestore(workspace=self.kwargs['wks'],
                                                       store_name=postgis_store_name,
                                                       pg_table=layer_name)
            layer_type = "vector"
        else:
            r = f'READ "no valid file extension: {file_extension} 400'
            layer_type = None
        return r, layer_type


_ = register_api(bp_geo, LayersAPI, "geo/layers", f"{app_api_base}/geo/layers/", pk="_id")
bp_geo.add_url_rule(app_api_base + '/geo/layers/<_id>.<string:_format>', view_func=_, methods=['GET'])


# view_func = LayerAPI.as_view("geo/layers")
# bp_geo.add_url_rule(f"{bcs_api_base}/geo/layers/", defaults={"_id": None}, view_func=view_func, methods=['GET', 'POST'])
# bp_geo.add_url_rule(f"{bcs_api_base}/geo/layers/<int:_id>", view_func=view_func, methods=['GET', 'PUT', 'DELETE'])


class RegionsAPI(MethodView):
    """
    Management of special "regions" layer through RESTful API
    """

    @n_session(read_only=True)
    def get(self, region_id=None):
        """
        ?? dede la GUI recibo el polígono marcado por el usuario
        @return:
        # dar lista de regiones para visualizar
        # entregar región a PDA
        """
        r = ResponseObject()
        db = g.n_session.db_session
        pg = g.n_session.postgis_db_session
        issues = []
        issues, geographic_region, count, status = get_content(db, GeographicRegion, issues, region_id)
        if status == 200 and geographic_region:
            if region_id is None:
                lines = False
                issues, regions, count, status = get_content(pg, Regions, issues, region_id)
            else:
                lines = True
                issues, regions, count, status = get_content(pg, Regions, issues, geographic_region.geo_id)

            bcs_df = pd.read_json(response_to_dataframe(geographic_region), lines=lines)
            postgis_df = pd.read_json(response_to_dataframe(regions), lines=lines)
            bcs_df = bcs_df.set_index(bcs_df["geo_id"]).drop(columns=["geo_id", "uuid"])
            postgis_df = postgis_df.set_index(postgis_df["id"]).drop(columns=["uuid", "id"])
            content = pd.concat([bcs_df, postgis_df], axis=1, join="inner")
        else:
            content = None
        return ResponseObject(issues=issues, status=status, content=content, content_type="application/json",
                              count=count).get_response()

    @n_session()
    def post(self):
        db = g.n_session.db_session
        pg = g.n_session.postgis_db_session
        t = request.json
        geographic_region_schema = getattr(GeographicRegion, "Schema")()
        regions_schema = getattr(Regions, "Schema")()
        regions_data = get_json_from_schema(Regions, t)
        geographic_region_data = get_json_from_schema(GeographicRegion, t)
        regions = regions_schema.load(regions_data, instance=Regions())
        pg.add(regions)
        pg.flush()
        geographic_region = geographic_region_schema.load(geographic_region_data, instance=GeographicRegion())
        geographic_region.uuid = regions.uuid
        geographic_region.geo_id = regions.id
        geographic_region.identity_id = g.n_session.identity.id
        db.add(geographic_region)
        db.flush()
        return ResponseObject(content=geographic_region, content_type="application/json").get_response()

    @n_session()
    def delete(self, region_id=None):
        issues = []
        db = g.n_session.db_session
        pg = g.n_session.postgis_db_session
        issues, geographic_region, count, status = get_content(db, GeographicRegion, issues, region_id)
        if status == 200 and geographic_region:
            issues, region, count, status = get_content(pg, Regions, issues, geographic_region.geo_id)
            db.delete(geographic_region)
            pg.delete(region)
        return ResponseObject(issues=issues, status=status).get_response()

    @n_session()
    def put(self, region_id=None):
        issues = []
        db = g.n_session.db_session
        pg = g.n_session.postgis_db_session
        t = request.json
        geographic_region_schema = getattr(GeographicRegion, "Schema")()
        regions_schema = getattr(Regions, "Schema")()
        issues, geographic_region, count, status = get_content(db, GeographicRegion, issues, region_id)
        if status == 200 and geographic_region:
            issues, region, count, status = get_content(pg, Regions, issues, geographic_region.geo_id)
            geographic_region = geographic_region_schema.load(get_json_from_schema(GeographicRegion, t),
                                                              instance=geographic_region)
            regions = regions_schema.load(get_json_from_schema(Regions, t), instance=region)
            db.add(geographic_region)
            pg.add(regions)
        return ResponseObject(issues=issues, status=status, count=count).get_response()


register_api(bp_geo, RegionsAPI, "geo/regions", f"{app_api_base}/geo/regions/", pk="region_id")


class StylesAPI(MethodView):
    """
    Create Styles:

    # 1. Dynamic styles from rasters coverages specifying color_Ramp
        input: color_ramp : dictionary,  matplotlib colormaps
        (https://matplotlib.org/3.3.0/tutorials/colors/colormaps.html), or list of colors
        c_ramp_1 = {
            'label 1 value': '#ffff55',
            'label 2 value': '#505050',
            'label 3 value': '#404040',
            'label 4 value': '#333333'
        }

        c_ramp_2 = 'RdYlGn'
        c_ramp_3 = [#ffffff, #453422,  #f0f0f0, #aaaaaa]

        geo.create_coveragestyle(raster_path=r'path\to\raster\file.tiff',
                            style_name='style_2',
                            workspace='demo',
                            color_ramp=c_ramp_,
                            cmap_type='values')

    # 2. Feature Styles:
        - Outline featurestyle: change boundary
        geo.create_outline_featurestyle(style_name='new_style' color="#3579b1" geom_type='polygon', workspace='demo')
        - Categorized featurestyle:
        geo.create_categorized_featurestyle(style_name='
        ', column_name='name_of_column', column_distinct_values=[1,2,3,4,5,6,7], workspace='demo')
        - Classified feature style:
        geo.create_classified_featurestyle(style_name='name_of_style' column_name='name_of_column', column_distinct_values=[1,2,3,4,5,6,7], workspace='demo')
    """

    @n_session()
    def get(self):
        """
        apply an style in a layer and publish
        @param id:
        @return:
        """
        # geo.publish_style(layer_name=layer, style_name='sld_file_name', workspace=style_name,
        #                 sld_version='1.0.0')# version?
        from biobarcoding.geo import geoserver_session
        styles = geoserver_session.get_styles()

    @n_session()
    def put(self, rampa):
        """
        create a new stylein geoserver
        @param rampa:
        @return:
        """
        # create new style
        # geo.upload_style(path=r'path\to\sld\file.sld', workspace='demo')
        pass

    @n_session()
    def post(self):
        pass

    def post_style_from_raster_file(self):
        from biobarcoding.geo import geoserver_session
        c_ramp = {
            'label 1 value': '#ffff55',
            'label 2 value': '#505050',
            'label 3 value': '#404040',
            'label 4 value': '#333333'
        }
        geoserver_session.create_coveragestyle(raster_path=r'path\to\raster\file.tiff',
                                               style_name='style_2',
                                               workspace='demo',
                                               color_ramp=c_ramp,
                                               cmap_type='values')

    def post_raster_style(self):
        pass

    def post_vector_style(self):
        from biobarcoding.geo import geoserver_session
        # - Outline featurestyle: change boundary
        geoserver_session.create_outline_featurestyle(style_name='new_style', color="#3579b1", geom_type='multipolygon',
                                                      workspace='ngd')
        # - Catagorized featurestyle:
        geoserver_session.create_catagorized_featurestyle(style_name="new",
                                                          column_name='name_of_column',
                                                          column_distinct_values=[1, 2, 3, 4, 5, 6, 7],
                                                          workspace='ngd',
                                                          color_ramp="tab20")
        # - Classified featurestyle:
        geoserver_session.create_classified_featurestyle(style_name='name_of_style', column_name='name_of_column',
                                                         column_distinct_values=[1, 2, 3, 4, 5, 6, 7], workspace='ngd')
