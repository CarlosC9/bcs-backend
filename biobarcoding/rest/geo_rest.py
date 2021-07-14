import io

import geopandas
import marshmallow.exceptions
import pandas as pd
from flask import Blueprint, request, g
from flask.views import MethodView
from biobarcoding.authentication import bcs_session
from biobarcoding.geo import geoserver_session, workspace_names, postgis_store_name
from biobarcoding.rest import bcs_api_base, ResponseObject, Issue, IType, register_api
from biobarcoding.db_models.geographics import GeographicRegion, Regions, GeographicLayer
import geopandas as gpd
import json
import regex as re
import requests
import pathlib
import numpy as np
import tempfile

from biobarcoding.rest.file_manager import FilesAPI

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


def geoserver_response(response) -> Issue:
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


def get_content(session, feature_class, issues, id_=None):
    content = None
    count = 0
    try:
        content = session.query(feature_class)
        if id_:
            content = content.filter(feature_class.id == id_).first()
        else:
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


# TODO Seems specific to BIOTA
def import_file_as_geojson(path):
    def f(x):
        return x[["RIQUEZA", "RAREZALOCA", "RAREZAINSU", "RAREZAREGI", "CODIGOTAX", "DENOMTAX"]].to_json(
            orient="records")

    gdf = gpd.read_file(path)
    for c in ["IDCELDA", "CODIGOTAX"]:
        gdf[c] = gdf[c].astype(np.int64)
    tmp = dict(zip(gdf["IDCELDA"], gdf["geometry"]))
    df = gdf.groupby("IDCELDA").apply(f).to_frame("taxa").reset_index()
    new_gdf = gpd.GeoDataFrame(df, crs=gdf.crs, geometry=df["IDCELDA"].replace(tmp))
    return new_gdf


class RegionsAPI(MethodView):
    """
    Management of special "regions" layer through RESTful API
    """

    @bcs_session(read_only=True)
    def get(self, region_id=None):
        """
        ?? dede la GUI recibo el polígono marcado por el usuario
        @return:
        # dar lista de regiones para visualizar
        # entregar región a PDA
        """
        r = ResponseObject()
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
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

    @bcs_session()
    def post(self):
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
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
        geographic_region.identity_id = g.bcs_session.identity.id
        db.add(geographic_region)
        db.flush()
        return ResponseObject(content=geographic_region, content_type="application/json").get_response()

    @bcs_session()
    def delete(self, region_id=None):
        issues = []
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
        issues, geographic_region, count, status = get_content(db, GeographicRegion, issues, region_id)
        if status == 200 and geographic_region:
            issues, region, count, status = get_content(pg, Regions, issues, geographic_region.geo_id)
            db.delete(geographic_region)
            pg.delete(region)
        return ResponseObject(issues=issues, status=status).get_response()

    @bcs_session()
    def put(self, region_id=None):
        issues = []
        db = g.bcs_session.db_session
        pg = g.bcs_session.postgis_db_session
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


register_api(bp_geo, RegionsAPI, "geo/regions", f"{bcs_api_base}/geo/regions/", pk="region_id")


class LayersAPI(MethodView):
    """
    GET:    the list of layers info from Geographiclayer table on bcs
            publish a result of filtering a layer (with geoserver_session.publish_featurestore_sqlview(sql_string))

    POST:   Three post cases (at least):
            1. post all new vector layer from geoJSON format (PDA case)
            2. Import a Raster o Vector layer from a file (via command)
            3. Save a temporal layer created by geoserver_session.publish_featurestore_sqlview(sql_string)
            create and store in GeographicLayer table as item
            then in postgis as table (format)
            and then publish  in geoserver as new layer

    PUT:    changing name, attributes or geometry

    DELETE: delete completely from GeographicLayer table (bcs) postgis and geoserver
    """

    # some default
    kwargs = {}
    issues = []
    status = int()

    @bcs_session()
    def get(self, _id=None):
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
        curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt "$API_BASE_URL/geo/layers/"

        @param _id: ID of a layer; empty for ALL layers (if no query params specified)
        @return:
        """
        self.issues = []
        from biobarcoding.geo import geoserver_session
        layers = None
        _filter = request.args.get("filter")
        key_col = request.args.get("key_col")
        db = g.bcs_session.db_session
        if _id:  # A layer
            self.issues, layers, count, self.status = get_content(db, GeographicLayer, self.issues, _id)
            if layers and layers.is_deleted:
                layers = None
                _, self.status = self.issues.append(Issue(IType.INFO, f'no data available')), 200
        elif not _filter and not key_col:  # All layers
            self.issues, layers, count, self.status = get_content(db, GeographicLayer, self.issues)
            if layers:
                layers = list(filter(lambda x: (x.is_deleted is False), layers))
        elif _filter and key_col:  # Temporary layer
            # TODO Ensure Temporary layer is also registered as BCS GeographicLayer
            tmp_view = f'tmpview_{g.bcs_session.identity.id}'
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
            layers = dict(tmpview=tmp_view, sql=sql, key_col=key_col)
        else:  # No information to elaborate a response
            _, self.status = self.issues.append(Issue(IType.ERROR, f"missing data")), 400

        return ResponseObject(issues=self.issues, status=self.status, content=layers).get_response()

    @bcs_session()
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
        """
        db = g.bcs_session.db_session
        self.issues = []
        system_layer = True
        self.kwargs = dict(name="unnamed")
        self._update_kwargs_from_request_data()
        self.kwargs["wks"] = workspace_names[0] if system_layer else workspace_names[1]
        # Put self.kwargs into a geographic layer
        geographic_layer_schema = getattr(GeographicLayer, "Schema")()
        geographic_layer_data = get_json_from_schema(GeographicLayer, self.kwargs)
        # Create object in memory
        geographic_layer = geographic_layer_schema.load(geographic_layer_data, instance=GeographicLayer())
        geographic_layer.identity_id = g.bcs_session.identity.id
        # Persist object in BCS DB
        db.add(geographic_layer)
        db.flush()
        # Receive file
        self._receive_and_prepare_file()
        # Store file in PostGIS
        layer_name = f"layer_{geographic_layer.id}"  # (internal) Geoserver layer name

        status = self._post_in_postgis(layer_name)
        # TODO Delete temporary file

        # Publish layer in Geoserver
        if status == 200:
            geographic_layer.in_postgis = True
        status, layer_type = self._publish_in_geoserver(layer_name)
        if status == 200:
            geographic_layer.published = True
            geographic_layer.geoserver_name = layer_name
        geographic_layer.layer_type = layer_type
        db.flush()
        return ResponseObject(issues=self.issues, status=self.status, content=geographic_layer).get_response()

    @bcs_session()
    def put(self, _id=None):
        self.issues = []
        self.kwargs = {}
        db = g.bcs_session.db_session
        # no puedo tener valores por defecto en el post
        self._update_kwargs_from_request_data()
        geographic_layer_schema = getattr(GeographicRegion, "Schema")()
        self.issues, geographic_layer, count, self.status = get_content(db, GeographicLayer, self.issues, _id)
        if geographic_layer:
            try:
                geographic_layer = geographic_layer_schema.load(self.kwargs, instance=geographic_layer)
            except marshmallow.exceptions.ValidationError as field:
                geographic_layer_data = get_json_from_schema(GeographicLayer, self.kwargs)
                geographic_layer = geographic_layer_schema.load(geographic_layer_data, instance=geographic_layer)
            if self.kwargs.get("data") or request.files:
                layer_name = f"layer_{geographic_layer.id}"
                geographic_layer.published = False
                geographic_layer.in_postgis = False
                db.flush()
                if request.files:
                    self._receive_and_prepare_file()
                if not self.kwargs.get("wks"):
                    self.kwargs["wks"] = geographic_layer.wks
                status = self._post_in_postgis(layer_name)
                if status == 200:
                    geographic_layer.in_postgis = True
                # re-publish in geoserver is need when the layer is changed (why?)
                # (note that it is not necessary when sql_view layer)
                status, layer_type = self._publish_in_geoserver(layer_name)
                if status == 200:
                    geographic_layer.published = True
                    geographic_layer.layer_type = layer_type
                    db.flush()
                else:
                    db.rollback()
        return ResponseObject(content=geographic_layer, issues=self.issues, status=self.status).get_response()

    @bcs_session()
    def delete(self, _id):
        """
        delete layer from PostGIS
        ¿and geoserver in case of raster layer?
        @param _id:
        @return:
        """
        from biobarcoding.geo import geoserver_session
        from biobarcoding import postgis_engine
        self.issues = []
        db = g.bcs_session.db_session
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
        @param filter: filter build by the user in GUI
        @return: sql query for postgis
        """
        # key_col = _filter["key"]
        # sql = _filter["sql"]
        # _filter = _filter.encode('utf-8')
        return _filter

    def _post_in_postgis(self, layer_name):
        """
        Create (store) feature layer in PostGIS

        :param layer_name:
        :return: None if error, status if success
        """
        from biobarcoding import postgis_engine
        from biobarcoding.geo import geoserver_session
        if self.kwargs["path"] != "":
            df = self._import_vector_file()
            if not isinstance(df, geopandas.GeoDataFrame):
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
                df.to_postgis(layer_name, postgis_engine, if_exists="replace")
        except ValueError as e:
            _, self.status = self.issues.append(Issue(IType.INFO, f"Layer named as {layer_name} error {e}")), 400
        else:
            self.kwargs["data"] = "dummy"  # "add geoserver layer" function looks for any value in "data" to create a layer that is stored in PostGIS
            _, self.status = self.issues.append(Issue(IType.INFO, f"layer stored in geo database")), 200
        return self.status

    def _publish_in_geoserver(self, layer_name):
        """
        Publish a layer -already in a registered store- in Geoserver

        :param layer_name:
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
            print(r)
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
            r, layer_type = self._import_raster_file(layer_name)
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

    def _update_kwargs_from_request_data(self):
        if request.json:
            self.kwargs.update(request.json)
        elif request.form:
            if request.form.get("metadata"):
                self.kwargs.update(json.loads(request.form["metadata"]))
        elif request.values:
            t = request.values.to_dict()
            for key, item in t.items():
                if isinstance(item, str):
                    try:
                        item = json.loads(item)
                        self.kwargs.update(item)
                    except json.decoder.JSONDecodeError:
                        self.kwargs.update(t)

    def _receive_and_prepare_file(self):
        import os
        formats = [".tif", ".gpkg", ".zip"]
        path = ""
        folder = tempfile.mkdtemp(prefix="bcs_")
        from werkzeug.utils import secure_filename

        # Download to a local file if self.kwargs.data contains a FilesAPI path
        if self.kwargs.get("data"):
            ctype, contents = FilesAPI.get_file_contents(g.bcs_session, self.kwargs.get("data"))
            file = f"submitted_file.zip"
            file_path = os.path.join(folder, secure_filename(file))
            with open(file_path, "wb") as f:
                f.write(contents)
        elif len(request.files) > 0:
            files = request.files.to_dict(flat=False)
            for _, value in files.items():
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

    def _import_vector_file(self):
        path = self.kwargs["path"]
        if self.kwargs.get("convert_to"):
            if self.kwargs["convert_to"] == "geojson":
                gdf = import_file_as_geojson(path)
                return gdf
            else:
                self.issues.append(Issue(IType.ERROR, "No valid data to convert"))
                self.status = 400
                return None
        else:
            try:
                gdf = gpd.read_file(path)
                return gdf
            except:
                return None

    def _import_raster_file(self, layer_name):
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


register_api(bp_geo, LayersAPI, "geo/layers", f"{bcs_api_base}/geo/layers/", pk="_id")

# view_func = LayerAPI.as_view("geo/layers")
# bp_geo.add_url_rule(f"{bcs_api_base}/geo/layers/", defaults={"_id": None}, view_func=view_func, methods=['GET', 'POST'])
# bp_geo.add_url_rule(f"{bcs_api_base}/geo/layers/<int:_id>", view_func=view_func, methods=['GET', 'PUT', 'DELETE'])


class StylesAPI(MethodView):
    """
    Create Styles:

    # 1. Dynamic styles from rasters corverages specifying color_Ramp
        input: color_ramp : dictionary,  maplorlib colormaps
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

    @bcs_session()
    def get(self):
        """
        apply an style in a layer and publish
        @param id:
        @return:
        """
        # geo.publish_style(layer_name=layer, style_name='sld_file_name', workspace=style_name,
        #                 sld_version='1.0.0')# version?
        styles = geoserver_session.get_styles()
        pass

    @bcs_session()
    def put(self, rampa):
        """
        create a new stylein geoserver
        @param rampa:
        @return:
        """
        # create new style
        # geo.upload_style(path=r'path\to\sld\file.sld', workspace='demo')
        pass

    @bcs_session()
    def post(self):
        pass

    def post_style_from_raster_file(self):
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
