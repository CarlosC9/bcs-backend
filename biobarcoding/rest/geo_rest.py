import io
import json
import math
import os
import pathlib
import tempfile
import traceback
import urllib
from typing import List, Tuple, Union, Dict
from urllib.parse import urlparse
from zipfile import ZipFile, ZIP_DEFLATED

import geopandas as gpd
import numpy as np
import pandas as pd
import regex as re
import requests
import seaborn as sns
from flask import Blueprint, request, g, Response
from flask.views import MethodView
from marshmallow import EXCLUDE
from matplotlib.colors import rgb2hex
from pandas.core.indexes.numeric import Float64Index
from sqlalchemy import Integer, and_
from sqlalchemy.orm.attributes import flag_modified

from .. import get_global_configuration_variable, app_acronym
from ..authentication import n_session
from ..common import generate_json
from ..common.decorators import Memoize
from ..common.helpers import read_yaml
from ..db_models.core import data_object_type_id, Dataset, PortInProcessInstance, CProcessInstance, CProcess, \
    CProcessPort
from ..db_models.geographics import GeographicRegion, Regions, GeographicLayer
from ..db_models.sysadmin import PermissionType
from ..geo import workspace_names, postgis_store_name
from . import app_api_base, ResponseObject, Issue, IType, register_api, app_proxy_base, get_content, auth_filter, \
    filter_parse, parse_request_params
from ..services.files import get_file_contents
from ..geo.biota import read_biota_file, generate_pda_species_file_from_layer, import_pda_result

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


def generate_sld_file(style_name, rule):
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


def generate_sld_category_rule(column_name, style_name, color, value, value_title, geom_type="Polygon"):
    return """
                <sld:Rule>
                    <sld:Name>{1}</sld:Name>
                    <sld:Title>{5}</sld:Title>
                    <ogc:Filter xmlns:ogc="http://www.opengis.net/ogc">
                        <ogc:PropertyIsEqualTo>
                            <ogc:PropertyName>{0}</ogc:PropertyName>
                            <ogc:Literal>{4}</ogc:Literal>
                        </ogc:PropertyIsEqualTo>
                    </ogc:Filter>
                    <sld:{7}Symbolizer>
                        <sld:Fill>
                          <sld:CssParameter name="fill">{2}</sld:CssParameter>
                        </sld:Fill>
                        <sld:Stroke>
                          <sld:CssParameter name="stroke">{3}</sld:CssParameter>
                          <sld:CssParameter name="stroke-width">{6}</sld:CssParameter>
                          <sld:CssParameter name="stroke-linejoin">mitre</sld:CssParameter>
                        </sld:Stroke>
                    </sld:{7}Symbolizer>
                </sld:Rule>
            """.format(
        column_name,  # 0
        style_name,  # 1
        color,  # 2
        "#000000",  # 3
        value,  # 4
        value_title,  # 5
        1 if geom_type == "Polygon" else 2,  # 6
        geom_type  # 7
    )


def generate_sld_continuous_rule(column_name, style_name, color, left, right, first, geom_type="Polygon"):
    return """
        <sld:Rule>
            <sld:Name>{1}</sld:Name>
            <sld:Title>{7}</sld:Title>
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
            <sld:{9}Symbolizer>
                <sld:Fill>
                  <sld:CssParameter name="fill">{2}</sld:CssParameter>
                </sld:Fill>
                <sld:Stroke>
                  <sld:CssParameter name="stroke">{3}</sld:CssParameter>
                  <sld:CssParameter name="stroke-width">{8}</sld:CssParameter>
                  <sld:CssParameter name="stroke-linejoin">mitre</sld:CssParameter>
                </sld:Stroke>
            </sld:{9}Symbolizer>
        </sld:Rule>
        """.format(
        column_name,  # 0
        style_name,  # 1
        color,  # 2
        "#000000",  # 3
        left - (1e-4 if first else 0),  # 4. First left value discounts a bit due to precision in floats
        right,  # 5
        "OrEqualTo" if first else "",  # 6
        f"{'[' if first else '('}{left:.1f}, {right:.1f}]",  # 7
        1 if geom_type == "Polygon" else 2,  # 8
        geom_type  # 9
        )


def generate_category_sld_file(
        style_name: str,
        column_name: str,
        categories: List[str],
        color_palette: str,
        geom_type: str = "Polygon"
):
    """
    Read color palette.
    - Match categories of palette with input categories
    - Assign remaining colors
    :param style_name:
    :param column_name:
    :param categories:
    :param color_palette:
    :param geom_type:
    :return:
    """
    palette_hex, values = get_fixed_style(color_palette)
    rule = ""
    intersect = set(categories) & set(values)
    difference = set(categories) - set(values)
    palette_indices = set(range(len(palette_hex)))
    cont_rules = 0
    # Categories matching those in the palette
    for i, category in enumerate(categories):
        if category in intersect:
            # Find index of category in palette
            index = values.index(category)
            palette_indices.remove(index)
            category_title = category
            rule += generate_sld_category_rule(column_name, style_name, palette_hex[index], category, category_title, geom_type)
            cont_rules += 1
    # Remaining categories
    palette_indices = sorted(list(palette_indices))
    use_remaining_colors = True
    for i, category in enumerate(categories):
        if category in difference:
            # Find index of category in palette
            if use_remaining_colors:
                index = palette_indices[i % len(palette_indices)]
            else:
                index = i % len(palette_hex)
            category_title = category
            rule += generate_sld_category_rule(column_name, style_name, palette_hex[index], category, category_title, geom_type)
            cont_rules += 1

    if cont_rules == 1:  # Force at least two rules (a bug in Geoserver?)
        rule += generate_sld_category_rule(column_name, style_name, palette_hex[index], category, category_title, geom_type)

    generate_sld_file(style_name, rule)


def generate_fixed_ramp_sld_file(
        style_name: str,
        column_name: str,
        color_ramp: str,
        geom_type: str = "Polygon"
):
    # Read RGB values and intervals from the name
    palette_hex, intervals = get_fixed_style(color_ramp)
    rule = ""
    for i, color in enumerate(palette_hex):
        rule += generate_sld_continuous_rule(column_name, style_name, color,
                                             intervals[i][0], intervals[i][1],
                                             i == 0,
                                             geom_type)

    generate_sld_file(style_name, rule)


def generate_dynamic_ramp_sld_file(
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
    interval = diff / n
    palette = sns.color_palette(color_ramp, int(n))
    palette_hex = [rgb2hex(i) for i in palette]
    # interval = N/4
    # color_values = [{value: color} for value, color in zip(values, palette_hex)]
    # print(color_values)
    rule = ""
    for i, color in enumerate(palette_hex):
        rule += generate_sld_continuous_rule(column_name, style_name, color,
                                             min_value + interval * i, min_value + interval * (i + 1),
                                             i == 0,
                                             geom_type)

    generate_sld_file(style_name, rule)


def create_and_publish_style(gs_session,
                             wkspc,
                             layer_name,
                             attribute,
                             min_value,
                             max_value,
                             categories: List[str],
                             style_name,
                             color_ramp: str = "RdYlGn_r",
                             cmap_type: str = "ramp",
                             number_of_classes: int = 5,
                             overwrite: bool = False,
                             geom_type: str = "polygon"):
    if get_styles()[color_ramp]:  # Dynamic
        generate_dynamic_ramp_sld_file(style_name=style_name,
                                       column_name=attribute,
                                       values=[min_value, max_value],
                                       n_intervals=number_of_classes,
                                       color_ramp=color_ramp,
                                       geom_type=geom_type)
    else:  # Static
        if categories:
            generate_category_sld_file(style_name, attribute, categories, color_ramp, geom_type)
        else:
            generate_fixed_ramp_sld_file(style_name, attribute, color_ramp, geom_type)
    gs_session.upload_style("style.sld", style_name, wkspc, sld_version="1.0.0", overwrite=overwrite)
    os.remove("style.sld")
    gs_session.publish_style(layer_name, style_name, wkspc)


def read_geolayer(layer_name):
    from .. import postgis_engine
    try:
        gdf = gpd.read_postgis(f"{layer_name}", postgis_engine, geom_col="geometry")
    except:
        try:
            gdf = gpd.read_postgis(f"{layer_name}", postgis_engine, geom_col="geom")
        except:
            gdf = pd.read_sql(f"select * from {layer_name}", postgis_engine)
    return gdf


def export_geolayer(db_sess, layer_id: int, layer_name: str, format_: str) -> Tuple[object, str]:
    from .. import postgis_engine
    if format_.lower() in ["geojson", "shp", "gpkg"]:
        try:
            gdf = gpd.read_postgis(f"{layer_name}", postgis_engine, geom_col="geometry")
        except:
            gdf = gpd.read_postgis(f"{layer_name}", postgis_engine, geom_col="geom")
    else:
        gdf = pd.read_sql(f"select * from {layer_name}", postgis_engine)
        # REMOVE GEOM column if it exists
        found = False
        for col in ["geom", "geometry"]+gdf.columns.to_list():
            if col in gdf.columns:
                try:
                    _ = _df_to_geodf(gdf, geom_col=col)
                    found = True
                    del gdf[col]
                    break
                except:
                    pass

    # Write to temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        if format_.lower() == "geojson":
            fname = f"{tmpdirname}{os.sep}{layer_name}.geojson"
            gdf.to_file(fname, driver="GeoJSON")
            content_type = "application/geo+json"
        elif format_.lower() == "shp":
            fname = f"{tmpdirname}{os.sep}{layer_name}.shp"
            gdf.to_file(fname)

            def zip_dir(zip_name: str, source_dir: Union[str, os.PathLike]):
                src_path = pathlib.Path(source_dir).expanduser().resolve(strict=True)
                with ZipFile(zip_name, 'w', ZIP_DEFLATED) as zf:
                    for file in src_path.rglob('*'):
                        if not str(file).endswith(".zip"):
                            zf.write(file, file.relative_to(src_path))

            fname = f"{tmpdirname}{os.sep}{layer_name}.zip"
            zip_dir(fname, tmpdirname)
            content_type = "application/zip"
        elif format_.lower() == "gpkg":
            fname = f"{tmpdirname}{os.sep}{layer_name}.gpkg"
            gdf.to_file(fname, driver="GPKG")
            content_type = "application/geopackage+sqlite3"
        elif format_.lower() == "csv":
            fname = f"{tmpdirname}{os.sep}{layer_name}.csv"
            gdf.to_csv(fname, index=False)
            content_type = "text/csv"
        elif format_.lower() == "xlsx":
            fname = f"{tmpdirname}{os.sep}{layer_name}.xlsx"
            gdf.to_excel(fname, index=False, sheet_name=layer_name, engine='xlsxwriter')
            content_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        # Read into memory
        with open(fname, 'rb') as fh:
            buf = io.BytesIO(fh.read())

    return buf.getvalue(), content_type


def _df_to_geodf(df, geom_col="geom", crs=None):
    """ COPIED FROM MODULE geopandas.io.sql """
    import shapely
    from geopandas import GeoDataFrame

    if geom_col not in df:
        raise ValueError("Query missing geometry column '{}'".format(geom_col))

    geoms = df[geom_col].dropna()

    if not geoms.empty:
        load_geom_bytes = shapely.wkb.loads
        """Load from Python 3 binary."""

        def load_geom_buffer(x):
            """Load from Python 2 binary."""
            return shapely.wkb.loads(str(x))

        def load_geom_text(x):
            """Load from binary encoded as text."""
            return shapely.wkb.loads(str(x), hex=True)

        if isinstance(geoms.iat[0], bytes):
            load_geom = load_geom_bytes
        else:
            load_geom = load_geom_text

        df[geom_col] = geoms = geoms.apply(load_geom)
        if crs is None:
            srid = shapely.geos.lgeos.GEOSGetSRID(geoms.iat[0]._geom)
            # if no defined SRID in geodatabase, returns SRID of 0
            if srid != 0:
                crs = "epsg:{}".format(srid)

    return GeoDataFrame(df, crs=crs, geometry=geom_col)


def obtain_dataset_properties(dataset_name: str, gdf: gpd.GeoDataFrame, lc_attributes: bool):
    """
    Obtain the list of attributes of a layer
    :param dataset_name:
    :param gdf: DataFrame to analyze
    :param lc_attributes: True -> attribute names to lowercase
    :return:
    """

    impact_levels_semaphore = set(["1", "2", "3", "4", "5", "0"])
    _ = []
    for i, c in enumerate(gdf.columns):
        # lower(*id), lower(id*), lower(*codigo*), lower(coordx), lower(coordy), lower(geom), lower(geometry)
        lc = c.lower()
        # To lowercase column names if requested
        if lc_attributes:
            c = lc

        # Skip empty columns
        if sum(gdf[c].notna()) == 0:
            _.append(dict(name=c, type="empty"))
            continue

        # Obtain attribute metadata
        print(f"Analyzing column {c}")
        # METADATA (per attribute)
        p_type = None  # Int, Numeric, String, Time, Category, Geometry
        p_type_2 = "measure"  # Dimension, measure or attribute
        # For numeric, the range of values
        min_v = None
        max_v = None
        # For category, the values
        categories = None
        p_subtype = None  # Only for "category" type, "int" or "string"
        # What [unit]
        nature = None
        unit = None
        # For representable attributes, the style name and colormap
        style_name = None
        colormap = None
        # Obtain the type of the attribute from the data
        if gdf.dtypes[i] == np.int64 or gdf.dtypes[i].name == "category":
            # If a column has a relatively small number of unique values, it is considered as a category
            uniq = gdf[c].unique()
            if uniq.size <= 10:
                p_type = "category"
                p_subtype = "int" if gdf.dtypes[i] == np.int64 else "string"
                # Find the set of different values
                if gdf.dtypes[i] == np.int64:
                    categories = sorted(list([int(x) for x in uniq]))
                else:
                    categories = sorted(list(uniq))
                # CATEGORY values must be strings (the SLD generation assumes this)
                categories = [str(x) for x in categories]
                # Find appropriate static style or create a dynamic one
                if set(uniq).issubset(impact_levels_semaphore) or len(uniq) < 6:
                    colormap = "__semaforo_impactos"
                    style_name = f"{dataset_name}_{c}"
            else:
                if gdf.dtypes[i] == np.int64:
                    p_type = "int"
                    min_v = min(gdf[c])
                    max_v = max(gdf[c])
                else:
                    p_type = "string"
        elif gdf.dtypes[i] in (np.float, np.float32, np.float64):
            # Numeric column
            tmp = gdf[c].values
            try:
                min_v = float(np.nanmin(tmp))
                max_v = float(np.nanmax(tmp))
                p_type = "numeric"
                if not math.isnan(min_v) and not math.isnan(max_v):
                    style_name = f"{dataset_name}_{c}"
                    colormap = "RdYlGn_r"
            except:
                p_type = "string"
            # TODO Find style: dynamic or static
        elif str(gdf.dtypes[i]) == "geometry":
            p_type = "geometry"
        else:
            # String column, do not create style
            p_type = "string"

        # Append attribute
        d = None
        if p_type in ("numeric", "int"):
            d = dict(name=c, col_type=p_type_2, data_type=p_type, min=min_v, max=max_v)
            if c not in ("id",):
                if style_name:
                    d["style"] = style_name
                if colormap:
                    d["colormap"] = colormap
        elif p_type == "category":
            d = dict(name=c, col_type=p_type_2, data_type=p_type, cat_type=p_subtype, categories=categories)
            if style_name:
                d["style"] = style_name
            else:
                d = dict(name=c, col_type=p_type_2, data_type=p_type)
            if colormap:
                d["colormap"] = colormap
        else:
            d = dict(name=c, col_type=p_type_2, data_type=p_type)
        if d:
            d["representable"] = "style" in d
            d["nature"] = nature
            d["unit"] = unit
            _.append(d)

    return _


def create_properties_and_geoserver_styles(properties: List[Dict],
                                           wks: str,
                                           layer_name: str,
                                           create_style=True,
                                           geometry_type: str = "Polygon"):
    """

    :param gdf: GeoDataFrame to analyze
    :param wks: Workspace that will container the style
    :param layer_name: Name of the layer, for the style
    :param lc_attributes: Lowercase attributes
    :param create_style: True to create a style, False to skip style creation (for non-geometry datasets)
    :return:
    """
    if create_style:
        from ..geo import geoserver_session
        for d in properties:
            style_name = d.get("style")
            if style_name:
                create_and_publish_style(geoserver_session,
                                         wkspc=wks,
                                         layer_name=layer_name,
                                         attribute=d["name"],
                                         min_value=d.get("min"), max_value=d.get("max"),
                                         number_of_classes=7,
                                         categories=d.get("categories"),
                                         color_ramp=d.get("colormap"),
                                         style_name=style_name,
                                         geom_type=geometry_type)


def refresh_properties_if_needed(db_sess, layer):
    # Analyze properties and if they are outdated, refresh
    refresh = False
    for d in layer.properties:
        if "representable" not in d:
            refresh = True
            break
    if refresh:
        gdf = read_geolayer(layer.geoserver_name)
        layer.properties = obtain_dataset_properties(layer.name, gdf, lc_attributes=True)
        db_sess.add(layer)
        db_sess.commit()


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

    def _export(self, sess, lay: GeographicLayer, _format) -> Tuple[object, str]:
        """

        :param sess:
        :param lay:
        :param _format:
        :return: A tuple with the exported layer and the content type
        """
        supported_formats = {"gpkg": ("Geopackage", export_geolayer),
                             "shp": ("Shapefile (zipped)", export_geolayer),
                             "geojson": ("GeoJSON", export_geolayer),
                             "csv": ("CSV (layers with Geometry column have it removed)", export_geolayer),
                             "xlsx": ("XLSX (layers with Geometry column have it removed)", export_geolayer),
                             "nexus": ("Nexus for PDA", generate_pda_species_file_from_layer),
                             "pda_simple": ("PDA simple", generate_pda_species_file_from_layer),
                             "species": ("List of species", generate_pda_species_file_from_layer),
                             "species_canon": ("List of normalized species", generate_pda_species_file_from_layer),
                             "pda_species": ("List of normalized and cleaned species names", generate_pda_species_file_from_layer)}
        if _format not in supported_formats.keys():
            self.issues.append(
                Issue(IType.ERROR, f'Could not export layer {lay.name} (internal name "{lay.geoserver_name}") '
                                   f'to format "{_format}". Supported: {", ".join(supported_formats.keys())}',
                      f"Export {lay.name}"))
            return None
        _, content_type = supported_formats[_format][1](sess, lay.id, lay.geoserver_name, _format)
        if _ is None:
            self.issues.append(
                Issue(IType.ERROR, f'Could not export layer {lay.name} (internal name "{lay.geoserver_name}").',
                      f"Export {lay.name} as {supported_formats[_format][0]}"))
            return None, None
        elif _ == "":
            self.issues.append(Issue(IType.ERROR,
                                     f'Empty layer {lay.name} (internal name "{lay.geoserver_name}").',
                                     f"Export {lay.name}"))
            return None, None
        else:
            return _, content_type

    @n_session(read_only=True)
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
        curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/"
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/?filter=%7B%22tags%22%3A%20%22biota%22%7D"

        (Metadata)
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/1"
        (Contents)
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/1.gpkg"
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/1.shp"
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/1.geojson"
        (Contents, only layers without Geometry column:)
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/1.csv"
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/1.xlsx"

        Specially, a Biota layer can be exported as one of 4 supported formats. The first two are for PDA, 3rd and 4th
        just enumerate available species:

        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/1.pda_simple"
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/1.nexus"
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/1.species"
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/1.species_canon"

        @param _id: ID of a layer; empty for ALL layers (if no query params specified)
        @param _format: export format
        @return:
        """

        def custom_geolayers_filter(filter, session=None):
            """
            for_geoprocesses, for_port_types

            sess.query(PortType).
            :param filter:
            :return:
            """
            clauses = []
            close_connection = session.transaction is None
            conn = session.connection()
            ids = None
            avoid_no_entity = False

            json_attributes = []
            for q_field, p_field, op in [("subjects", "themes_id", "?|"),
                                         ("sources", "source_id", "="),
                                         ("crs", "crs_id", "=")]:
                if filter.get(q_field):
                    # Thanks to (?| does not support a list of integers, only of strings, so a workaround):
                    # https://stackoverflow.com/questions/48321403/postgres-jsonb-int-array-contains-any-of-array-values
                    if op == "?|":
                        _ = []
                        for m in filter[q_field]["unary"]:
                            _.append(f"ds.attributes->'{p_field}' @> '{int(m)}'::jsonb")
                        json_attributes.append(f'({" OR ".join(_)})')
                    else:
                        json_attributes.append(
                            f"ds.attributes->'{p_field}' {op} '{int(filter[q_field]['unary'])}'::jsonb")
            if json_attributes:
                sql_select = f"SELECT DISTINCT(ds.id) " \
                             f"FROM {app_acronym}_fos ds " \
                             f"WHERE {' AND '.join(json_attributes)}"
                sql_result = conn.execute(sql_select)
                _ = [r[0] for r in sql_result]
                if not _:
                    avoid_no_entity = True
                if ids:
                    ids.intersect(_)
                else:
                    ids = set(_)

            for q_field, port_direction_condition in [("used_as_input_of", "gp.input"),
                                                      ("resulting_from", "not gp.input")]:
                if filter.get(q_field):
                    _ = ', '.join(['\'' + str(i) + '\'' for i in ids]) if ids else ''
                    remaining_layers_condition = f" AND ds.id IN ({_})" if ids else ""
                    lst = ', '.join([f"{ds}" for ds in filter[q_field]["unary"]])
                    sql_select = f"SELECT DISTINCT(ds.id) " \
                                 f"FROM {Dataset.__tablename__} ds " \
                                 f"JOIN {PortInProcessInstance.__tablename__} ppi ON ds.id=ppi.dataset_id " \
                                 f"JOIN {CProcessInstance.__tablename__} pi ON ppi.process_instance_id=pi.id " \
                                 f"JOIN {CProcess.__tablename__} g ON pi.instantiated_process_id=g.id " \
                                 f"JOIN {CProcessPort.__tablename__} gp ON gp.id=ppi.port_id " \
                                 f"WHERE {port_direction_condition} AND g.id IN ({lst})" \
                                 f"{remaining_layers_condition}"
                    sql_result = conn.execute(sql_select)
                    _ = [r[0] for r in sql_result]
                    if not _:
                        avoid_no_entity = True
                    if ids:
                        ids.intersect(_)
                    else:
                        ids = set(_)

            if filter.get('for_geoprocesses'):  # May be used as input of a geoprocess (it could have been used already)
                _ = ', '.join(['\'' + str(i) + '\'' for i in ids]) if ids else ''
                remaining_layers_condition = f" AND ds.id IN ({_})" if ids else ""
                lst = ', '.join([f"{ds}" for ds in filter["for_geoprocesses"]["unary"]])
                sql_select = f"SELECT DISTINCT(ds.id) " \
                             f"FROM {app_acronym}_datasets ds JOIN {app_acronym}_dataset_port_types ds_pt ON ds.id=ds_pt.dataset_id " \
                             f"WHERE ds_pt.port_type_id IN " \
                             f"(SELECT DISTINCT(gp.port_type_id) " \
                             f" FROM {app_acronym}_processes g JOIN {app_acronym}_processes_ports gp ON g.id=gp.process_id " \
                             f" WHERE gp.input AND g.id IN ({lst})) " \
                             f"{remaining_layers_condition}"
                sql_result = conn.execute(sql_select)
                _ = [r[0] for r in sql_result]
                if not _:
                    avoid_no_entity = True
                if ids:
                    ids.intersect(_)
                else:
                    ids = set(_)

            if filter.get('for_port_types'):
                _ = ', '.join(['\'' + str(i) + '\'' for i in ids]) if ids else ''
                remaining_layers_condition = f" AND ds.id IN ({_})" if ids else ""
                lst = ', '.join([f"{ds}" for ds in filter["for_port_types"]["unary"]])
                sql_select = f"SELECT DISTINCT(ds.id) " \
                             f"FROM {app_acronym}_datasets ds JOIN {app_acronym}_dataset_port_types ds_pt ON ds.id=ds_pt.dataset_id " \
                             f"WHERE ds_pt.port_type_id IN ({lst}) " \
                             f"{remaining_layers_condition}"
                sql_result = conn.execute(sql_select)
                _ = [r[0] for r in sql_result]
                if not _:
                    avoid_no_entity = True
                if ids:
                    ids.intersect(_)
                else:
                    ids = set(_)

            if filter.get("case_studies_"):
                _ = ', '.join(['\'' + str(i) + '\'' for i in ids]) if ids else ''
                remaining_layers_condition = f" AND ds.id IN ({_})" if ids else ""
                lst = ', '.join([f"{ds}" for ds in filter["case_studies_"]["unary"]])
                sql_select = f"SELECT DISTINCT(ds.id) " \
                             f"FROM {app_acronym}_datasets ds " \
                             f"JOIN {app_acronym}_case_studies_functional_objects cs ON ds.id=cs.functional_object_id " \
                             f"WHERE cs.case_study_id IN ({lst}) " \
                             f"{remaining_layers_condition}"
                sql_result = conn.execute(sql_select)
                _ = [r[0] for r in sql_result]
                if not _:
                    avoid_no_entity = True
                if ids:
                    ids.intersect(_)
                else:
                    ids = set(_)

            if avoid_no_entity and len(ids) == 0:
                ids = [-1]

            if ids and len(ids) > 0:
                clauses.append(GeographicLayer.id.in_(list(ids)))

            if close_connection:
                conn.close()

            return clauses

        # ---------------------------------------------------------------------
        self.issues = []
        layer = None
        count = 1
        key_col = request.args.get("key_col")
        db_sess = g.n_session.db_session
        if _id:  # One layer or Job (return layer associated to job)
            if _id.startswith("job"):  # Layers associated with a job
                job_id = _id[len("job"):]
                # TODO Return more than one layer if there are more than one
                layer = db_sess.query(GeographicLayer). \
                    filter(GeographicLayer.attributes["job_id"].astext.cast(Integer) == int(job_id)).first()
                if layer:
                    self.status = 200
                else:
                    self.status = 400
            else:
                self.issues, layer, count, self.status = get_content(db_sess, GeographicLayer, self.issues, _id)
                refresh_properties_if_needed(db_sess, layer)  # Read layer from PostGIS and update properties
            if layer and layer.is_deleted:  # Soft deleted
                layer = None
                _, self.status = self.issues.append(Issue(IType.INFO, f'no data available')), 200
            else:
                if _format:  # Download Layer contents
                    content, content_type = self._export(db_sess, layer, _format=_format)
                    if content:  # Export layer
                        return Response(content, mimetype=content_type, status=200)
                else:  # Layer metadata
                    if layer:
                        serializer = layer.Schema()
                        serializer.dump(layer)
        elif not key_col:  # ALL LAYERS layers (maybe filtered)
            # ACL Filter, here !!!
            purpose_id = db_sess.query(PermissionType).filter(PermissionType.name == "read").one().id
            ids_clause = db_sess.query(GeographicLayer.id). \
                filter(auth_filter(GeographicLayer, purpose_id, [data_object_type_id['geolayer']],
                                   identity_id=g.n_session.identity.id)).subquery()
            query = db_sess.query(GeographicLayer).filter(GeographicLayer.id.in_(ids_clause))
            # Modify request.values to filter by "is_deleted" attribute of GeographicLayer
            _ = dict(request.values)
            if "filter" in _:
                f = json.loads(urllib.parse.unquote(_["filter"]))
            else:
                f = {}
            if isinstance(f, dict):
                f["is_deleted"] = {"op": "eq", "unary": "False"}
            _["filter"] = urllib.parse.quote(json.dumps(f))
            request.values = _

            # Get content
            self.issues, layer, count, self.status = get_content(db_sess, GeographicLayer, self.issues,
                                                                 aux_filter=custom_geolayers_filter,
                                                                 query=query)
            if layer:
                # Modify "wms_url" for local layers
                tmp = urlparse(request.base_url)
                # WORKAROUND - tmp.scheme is always "http" in spite the request being https. So assume "https"
                if not tmp.netloc.lower().startswith(("localhost", "172.17.0.1", "127.0.0.1", "0.0.0.0", "10.141")):
                    scheme = "https"
                else:
                    scheme = tmp.scheme
                base_url = f"{scheme}://{tmp.netloc}"
                # Rewritten to avoid modification of geographic layers that would generate ORM events
                enh_content = []
                for inst in layer:
                    _ = json.loads(generate_json(inst))
                    if _.get("wms_url", "").endswith(f"{app_proxy_base}/geoserver/wms"):
                        _["wms_url"] = f"{base_url}{app_proxy_base}/geoserver/wms"
                    enh_content.append(_)
                layer = enh_content

        # elif _filter and key_col:  # Temporary layer
        #     # TODO Ensure Temporary layer is also registered as BCS GeographicLayer
        #     tmp_view = f'tmpview_{g.n_session.identity.id}'
        #     sql = self._create_sql(_filter)
        #     r = geoserver_session.delete_layer(layer_name=tmp_view, workspace=workspace_names[1])
        #     print(r)
        #     r = geoserver_session.publish_featurestore_sqlview(store_name=postgis_store_name, name=tmp_view, sql=sql,
        #                                                        key_column=key_col,
        #                                                        workspace=workspace_names[1])
        #     print(r)
        #     issue, self.status = geoserver_response(r)
        #     if not self._check_layer():
        #         _, self.status = self.issues.append(Issue(IType.ERROR, f"Error executing request for geoserver")), 500
        #     self.issues.append(issue)
        #     layer = dict(tmpview=tmp_view, sql=sql, key_col=key_col)
        else:  # No information to elaborate a response
            _, self.status = self.issues.append(Issue(IType.ERROR, f"missing data")), 400

        return ResponseObject(issues=self.issues, status=self.status, content=layer, count=count).get_response()

    @staticmethod
    def _exclude(c):
        """
        Check if column "c" has to be excluded from the list of attributes of a layer showed to users

        :param c:
        :return:
        """
        s = c.lower()
        return s.startswith("id") or s.endswith("id") or "codigo" in s or s in ["coordx", "coordy", "geom", "geometry"]

    def fill_types_and_case_studies(self, session, geographic_layer, geographic_layer_data):
        if "types" in geographic_layer_data:
            # Update new
            types = self.kwargs["types"]
            for i, t in enumerate(geographic_layer.types):
                t.dataset = geographic_layer
                t.port_type_id = types[i]
                session.add(t)
        if "case_studies" in geographic_layer_data:
            # Update new
            case_studies = self.kwargs["case_studies"]
            for i, t in enumerate(geographic_layer.case_studies):
                t.functional_object = geographic_layer
                t.case_study_id = case_studies[i]
                session.add(t)

    # Create the SQL from the request
    @staticmethod
    def _create_sql_from_request(session, query_json):
        """
        { "select": ["c1", "c2"],
          "from": ["t1", "t2],
          "join": [{"left": "", "right": "", "operator": ""}],
          "where": [],
          }

        :param query_json:
        :return:
        """
        def get_table_names(layer_codes):
            from ..db_models.geographics import GeographicLayer
            codes = [int(code) for code in layer_codes]
            layers = session.query(GeographicLayer).filter(GeographicLayer.id.in_(codes)).all()
            # TODO Knowing the layers and fields involved, it is possible to know the metadata of the original fields,
            #  without reading. This is the place where this information can be obtained.
            # fields = {}
            # for layer in layers:
            #     attributes = layer.attributes

            return {l.id: l.geoserver_name for l in layers}

        def get_table_field(fld):
            _ = fld.split("$")
            table = layer_names_dict[int(_[0])]
            field = _[1]
            return table, field, f"{table}.{field}"

        def select_clause(select_):
            select_a = []
            for r in select_["rules"].values():
                tab, fld, fld_sql = get_table_field(r["field"])
                if r["operator"] == "not_function":
                    _ = f"{fld_sql}"
                else:
                    params = ', '.join(r["value"])
                    _ = f"{r['operator']}({fld_sql}{params})"
                if r["alias"]:
                    _ += f" AS {r['alias']}"
                select_a.append(_)
            else:
                select_a.append("*")

            return ", ".join(select_a)

        def from_join_clause(t, join_):
            """
{
   "condition":"JOIN",
   "rules":[
      {
         "id":"7$x_coord",
         "field":"7$x_coord",
         "type":"string",
         "input":"select",
         "operator":"INNER",
         "value":"7$y_coord"
      },
      {
         "id":"7$x_coord",
         "field":"7$x_coord",
         "type":"string",
         "input":"select",
         "operator":"OUTER",
         "value":"6$riqueza"
      },
      {
         "id":"6$rarezaregi",
         "field":"6$rarezaregi",
         "type":"string",
         "input":"select",
         "operator":"INNER",
         "value":"1294$rcp85_tasmax_2040"
      }
   ],
   "valid":true
}
            :param t:
            :param join_:
            :return:
            """
            if len(t) == 1:
                return t[0]

            already_joined = set()
            from_join_str = ""
            for j in join_["rules"]:
                tab_l, fld_l, fld_sql_l = get_table_field(j["field"])
                tab_r, fld_r, fld_sql_r = get_table_field(j["value"])
                if tab_l not in already_joined:
                    if j["operator"].lower() in ("inner", "outer", "left", "right"):
                        from_join_str += f"{tab_l} {j['operator']} JOIN {tab_r} ON {fld_sql_l} = {fld_sql_r}"
                    else:  # Geom
                        from_join_str += f"{tab_l} JOIN {tab_r} ON {j['operator']}({fld_sql_l}, {fld_sql_r})"
                else:
                    if j["operator"].lower() in ("inner", "outer", "left", "right"):
                        from_join_str += f"{j['operator']} JOIN {tab_r} ON {fld_sql_l} = {fld_sql_r}"
                    else:  # Geom
                        from_join_str += f"JOIN {tab_r} ON {j['operator']}({fld_sql_l}, {fld_sql_r})"

                already_joined.add(tab_l)
                from_join_str

            return from_join_str

        def where_clause(where_):
            params = set()
            where_str = ""
            # Boolean operators
            if where_["condition"].lower() in ("and", "or", "not"):
                subs = []
                for rule in where_["rules"]:
                    _, p2 = where_clause(rule)
                    subs.append(_)
                    params.update(p2)
                if where_["condition"] == "not":
                    where_str = f"NOT ({subs[0]})"
                else:
                    where_str += "(" + f") {where_['condition']} (".join(subs) + ")"
            elif where_["operator"].lower() in ("equals", "equals_ignore_case",
                                                "=", ">", "<", ">=", "<=", "!=",
                                                "like", "ilike"):
                if where_["value"].startswith("$"):
                    params.add(where_["value"][1:])
                operator = where_["operator"].lower()
                if operator in ("equals_ignore_case", "equals"):
                    operator = "="
                if operator in ("equals_ignore_case", ):
                    function = "LOWER("
                else:
                    function = "("
                tab, fld, fld_sql = get_table_field(where_["field"])
                value = f"'{where_['value']}'" if where_["type"] == "string" else where_["value"]
                where_str += f"{function}{fld_sql}) {operator} {value}"
            elif where_["operator"].lower() == "in":
                if where_["value"].startswith("$"):
                    params.add(where_["value"][1:])
                tab, fld, fld_sql = get_table_field(where_["field"])
                values = ','.join([f"'{v}'" if where_["type"] == "string" else v for v in where_["value"]])
                where_str += f"{fld_sql} IN ({values})"

            return where_str, params

        # ---------------------------------------------------------------------

        layer_names_dict = get_table_names(query_json["layers"])
        select = select_clause(layer_names_dict, query_json["select"])
        from_ = from_join_clause(layer_names_dict, query_json["layers"], query_json["join"])
        where = where_clause(layer_names_dict, query_json["where"])

        return f"SELECT {select} FROM {from_} WHERE {where}"

    def _create_postgis_view_and_read(self, sql, layer_name, lower_case_attributes, read=False):
        """

        :param sql:
        :param layer_name:
        :param lower_case_attributes:
        :param read:
        :return:
        """
        from .. import postgis_engine
        # Create a view (DROP {materialized} VIEW if exists first)
        # It can be materialized or not
        cached = True
        materialized = "MATERIALIZED" if cached else ""
        postgis_engine.execute(f"DROP {materialized} VIEW IF EXISTS {layer_name}")
        s = f'CREATE {materialized} VIEW "{layer_name}" AS sql'
        postgis_engine.execute(sql)
        self.issues.append(Issue(IType.INFO, f"layer created and stored in geo database"))
        self.status = 200
        if read:
            # Read the data
            df = pd.read_sql(f"SELECT * FROM {layer_name}", postgis_engine)
            # TODO Try to convert to GeoPandas
        else:
            df = None

        return self.status, df

    @n_session()
    def post(self):
        """
        (CURRENTLY, NO) 1. import a raster file (tif and twf file)
        2. import vector layer: from SHP file (import shp shx...), GeoPackage, GeoJSON, CSV, ...
        3. create a new layer from a temp_view
        4. create a new vector file from shape file and convert data to json
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
        curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        curl --cookie app-cookies.txt -F "layer_file=@/home/rnebot/GoogleDrive/AA_NEXTGENDEM/plantae_canarias/Plantas.zip;type=application/zip" -F "metadata={\"name\": \"capa_1\", \"wks\": \"ngd\", \"attributes\": {\"tags\": [\"tag1\", \"tag2\"]}};type=application/json" "$API_BASE_URL/geo/layers/"

        POST from FilesAPI:
        First, upload a file to FilesAPI:
        curl --cookie app-cookies.txt -H "Content-Type: application/zip" -XPUT --data-binary @/home/rnebot/GoogleDrive/AA_NEXTGENDEM/plantae_canarias/Plantas.zip "$API_BASE_URL/files/f1/f2/f3/plantas.zip.content"

        Then, do the POST, like
        curl --cookie app-cookies.txt -XPOST "$API_BASE_URL/geo/layers/?filesAPI=%2Ff1%2Ff2%2Ff3%2Fplantas.zip&job_id=6"
        """
        def create_geographic_layer_entity(session, kwargs):
            geographic_layer_schema = getattr(GeographicLayer, "Schema")()
            geographic_layer_data = get_json_from_schema(GeographicLayer, self.kwargs)
            geographic_layer_data["id"] = 0
            # Create object in memory
            layer_ = geographic_layer_schema.load(geographic_layer_data, instance=GeographicLayer())
            layer_.id = None
            layer_.identity_id = g.n_session.identity.id
            self.fill_types_and_case_studies(session, layer_, geographic_layer_data)
            return layer_

        session = g.n_session.db_session
        self.issues = []
        self.kwargs = self._build_args_from_request_data()
        system_layer = True  # or user layer
        self.kwargs["wks"] = workspace_names[0] if system_layer else workspace_names[1]

        # Create layer entity, and persist it
        geographic_layer = create_geographic_layer_entity(session, self.kwargs)
        session.add(geographic_layer)
        session.flush()

        # Set URL
        tmp = urlparse(request.base_url)
        base_url = f"{tmp.scheme}://{tmp.netloc}"
        geographic_layer.wms_url = f"{base_url}{app_proxy_base}/geoserver/wms"

        lower_case_attributes = True
        layer_name = f"layer_{geographic_layer.id}"  # (internal) Geoserver layer name

        # Receive file or create a view from an SQL Query
        # NOTE: Creating a view can take a long time. Consider
        if self._file_posted():
            # Store file locally
            self._receive_and_prepare_file()
            # Store file in PostGIS
            status, gdf, has_geom_column = self._read_and_store_into_postgis(layer_name, lower_case_attributes)
        else:  # An SQL to create a new layer
            sql, has_geom_column = self._create_sql_from_request(session, self.kwargs["layer_from_query"])
            status, _ = self._create_postgis_view_and_read(sql, layer_name, lower_case_attributes)

        # Geometry type
        if not geographic_layer.attributes:
            geographic_layer.attributes = {}

        if not gdf.empty:
            geographic_layer.attributes["geom_type"] = "Line" if gdf.geom_type[0] == "MultiLineString" else "Polygon"
        else:
            geographic_layer.attributes["geom_type"] = "Polygon"  # Not known yet...

        # Delete temporary file
        os.remove(self.kwargs["path"])
        if status == 200:
            geographic_layer.in_postgis = True

        if has_geom_column:  # Only representable layers can be shown with GeoServer
            # Publish layer in Geoserver
            status, layer_type = self._publish_in_geoserver(layer_name)
            if status == 200:
                geographic_layer.published = True
                geographic_layer.geoserver_name = layer_name
                if not gdf.empty:
                    geographic_layer.properties = obtain_dataset_properties(layer_name, gdf, lower_case_attributes)
                    create_properties_and_geoserver_styles(geographic_layer.properties,
                                                           self.kwargs["wks"],
                                                           layer_name,
                                                           geometry_type=geographic_layer.attributes["geom_type"])
        else:
            layer_type = "no_explicit_geometry_layer"
            geographic_layer.geoserver_name = layer_name
            if not gdf.empty:
                geographic_layer.properties = obtain_dataset_properties(layer_name, gdf, lower_case_attributes)
                create_properties_and_geoserver_styles(geographic_layer.properties,
                                                       self.kwargs["wks"],
                                                       layer_name,
                                                       create_style=False)

        geographic_layer.layer_type = layer_type
        session.flush()

        # Update PortInProcessInstance if metadata is available. IMPORT PACKAGES LOCALLY
        if "attributes" in self.kwargs and \
                "port_id" in self.kwargs["attributes"] and \
                "instance_id" in self.kwargs["attributes"]:
            from ..db_models.core import PortInProcessInstance
            port_id = self.kwargs["attributes"]["port_id"]
            instance_id = self.kwargs["attributes"]["instance_id"]
            _ = session.query(PortInProcessInstance). \
                filter(and_(PortInProcessInstance.port_id == int(port_id),
                            PortInProcessInstance.process_instance_id == int(instance_id))).one()
            _.dataset_id = geographic_layer.id

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
        curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        curl --cookie app-cookies.txt -X PUT -F "metadata={\"name\": \"capa_12\", \"wks\": \"ngd\", \"attributes\": {\"tags\": [\"tag1\", \"tag2\", \"tag3\"]}};type=application/json" "$API_BASE_URL/geo/layers/"

        :param _id:
        :return:
        """
        session = g.n_session.db_session
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
        self.issues, geographic_layer, count, self.status = get_content(session, GeographicLayer, self.issues, _id)
        if geographic_layer:
            # Update layer metadata
            if "types" in geographic_layer_data:
                # Remove previous types
                for t in geographic_layer.types:
                    session.delete(t)
            if "case_studies" in geographic_layer_data:
                # Remove previous case studies
                for t in geographic_layer.case_studies:
                    session.delete(t)
            geographic_layer = geographic_layer_schema.load(geographic_layer_data,
                                                            instance=geographic_layer,
                                                            partial=True, unknown=EXCLUDE)
            self.fill_types_and_case_studies(session, geographic_layer, geographic_layer_data)

            # NOTE: Update layer content disabled (it should be another layer). Previous code available in Git
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
        from .. import postgis_engine
        from ..geo import geoserver_session
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
            if geographic_layer.in_postgis is False and geographic_layer.published is False:
                db.delete(geographic_layer)
        return ResponseObject(content=None, issues=self.issues, status=self.status).get_response()

    def _prepare_df(self, df, lower_columns):
        """
        Prepare dataframe columns:
         - cast to numeric and to category, when possible
         - replace <null> geometry with Empty Geometry
         - if specified, lowercase columns
        :param df: Dataframe to process
        :param lower_columns: True to put all columns in lowercase
        :return:
        """
        has_geom_column = df.geometry.dropna().shape[0] > 0

        # Convert columns to number, when possible
        for i, c in enumerate(df.columns):
            if self._exclude(c):
                continue

            # Skip empty columns
            if sum(df[c].notna()) == 0:
                continue

            uniq = df[c].unique()
            is_category = len(uniq) <= 10  # TODO
            is_dimension = len(uniq) <= 100  # TODO

            # Check if column "c" can be converted to "category" or to "number"
            if df.dtypes[i] in (np.int64, np.float64, np.dtype("O")) and uniq.size <= 10:
                df.loc[df[c].isna(), c] = "-" if df.dtypes[i] == np.dtype("O") else 0  # Replace NA with a null value
                df[c] = df[c].astype('category')  # Convert to category
                if isinstance(df[c].cat.categories, Float64Index):
                    # Try to convert to int64 categories
                    try:
                        df[c].cat.categories = df[c].cat.categories.astype("int64")
                    except:
                        pass
                else:
                    if "" in uniq:
                        df[c] = df[c].cat.rename_categories({'': '-'})  # Rename empty string category
            elif df.dtypes[i] == np.dtype("O"):
                # Numeric?
                try:
                    df[c] = pd.to_numeric(df[c])  # Convert to numeric
                except:
                    pass

        # Missing geometries substituted by empty geometries
        from shapely.geometry import Polygon
        df.geometry[df.geometry.isna()] = Polygon([])

        # All columns to lower case
        if lower_columns:
            df.columns = [c.lower() for c in df.columns]

        return has_geom_column

    def _read_and_store_into_postgis(self, layer_name, lower_columns=True):
        """
        Create (store) feature layer in PostGIS
        NOTE: implicit parameters in self.kwargs

        :param layer_name:
        :param lower_columns: rename all columns to lower case
        :return: None if error, status if success
        """
        from .. import postgis_engine

        # (1/3) READ
        if self.kwargs["path"] != "":
            df = self._read_vector_file()
            if not isinstance(df, gpd.GeoDataFrame):
                return None
        elif "data" in self.kwargs.keys():
            df = gpd.GeoDataFrame.from_features(self.kwargs["data"]['features'])
        # elif "layer_name" in self.kwargs.keys():
        #     layer = geoserver_session.get_layer(layer_name=self.kwargs["layer_name"])
        #     if isinstance(layer, dict):
        #         layer = geoserver_session.get_layer(layer_name=self.kwargs["layer_name"])
        #         # si no hay tmpview va a lanzar error
        #         if layer["layer"].get("type") == 'VECTOR':
        #             # TODO OPTION 1: export data and execute a Python function
        #             # OPTION 2: execute a SQL view statement into PostGIS
        #             df = self._get_df_from_view()
        #     else:  # go to geoserver
        #         return None

        # (2/3) PREPARE
        has_geom_column = self._prepare_df(df, lower_columns=lower_columns)

        # (3/3) STORE
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
                if has_geom_column:
                    # Write to PostGIS
                    df.to_postgis(layer_name, postgis_engine, if_exists="replace")
                else:
                    del df[df.geometry.name]
                    df.to_sql(layer_name, postgis_engine, if_exists="replace", index=False)
        except ValueError as e:
            _, self.status = self.issues.append(Issue(IType.INFO, f"Layer named as {layer_name} error {e}")), 400
        else:
            # "add geoserver layer" function looks for any value in "data" to create a layer that is stored in PostGIS
            self.kwargs["data"] = "dummy"
            _, self.status = self.issues.append(Issue(IType.INFO, f"layer stored in geo database")), 200
        return self.status, df, has_geom_column

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
        from ..geo import geoserver_session
        layer_type = "vector"
        if "data" in self.kwargs.keys():
            r = geoserver_session.delete_layer(layer_name=layer_name, workspace=self.kwargs['wks'])
            print(r)
            r = geoserver_session.publish_featurestore(workspace=self.kwargs['wks'],
                                                       store_name=postgis_store_name,
                                                       pg_table=layer_name)
            if attribute:
                create_and_publish_style(geoserver_session,
                                         wkspc=self.kwargs['wks'],
                                         layer_name=layer_name,
                                         attribute=attribute,
                                         min_value=min_value, max_value=max_value,
                                         categories=None,
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
        # elif self.kwargs.get("filter"):
        #     sql = self.kwargs["filter"]
        #     key_col = self.kwargs["key_col"]
        #     r = geoserver_session.publish_featurestore_sqlview(store_name=postgis_store_name,
        #                                                        name=layer_name,
        #                                                        sql=sql,
        #                                                        key_column=key_col,
        #                                                        workspace=self.kwargs['wks'])
        #     layer_type = "sqlview"
        elif request.files:
            r, layer_type = self._read_raster_file(layer_name)
        else:
            r = "No  valid geographic data specified, 400"

        issue, self.status = geoserver_response(r)
        self.issues.append(issue)
        return self.status, layer_type

    # def _tmpview_sql(self, layer_name):
    #     from ..geo import geoserver_session
    #     layer = geoserver_session.get_layer(layer_name=layer_name)
    #     if not isinstance(layer, dict):
    #         return None, None
    #     response = requests.get(layer["layer"]["resource"]["href"],
    #                             auth=(geoserver_session.username, geoserver_session.password))
    #     if response.status_code == 200:
    #         tmp_view_json = response.text
    #         tmp_view_dict = json.loads(tmp_view_json)
    #         sql = tmp_view_dict["featureType"]["metadata"]["entry"]["virtualTable"]["sql"].rstrip("\n")
    #         key_col = tmp_view_dict["featureType"]["metadata"]["entry"]["virtualTable"]["keyColumn"]
    #         return sql, key_col
    #     else:
    #         return None, None
    #
    # def _get_df_from_view(self):
    #     sql, key_col = self._tmpview_sql(self.kwargs["layer_name"])
    #     if sql:
    #         from .. import postgis_engine
    #         df = gpd.GeoDataFrame.from_postgis(sql, postgis_engine)
    #         return df
    #     else:
    #         _, self.status = self.issues.append(Issue(IType.INFO, f"no tmpview available")), 400
    #         return None

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
            else:
                args["name"] = "Nombre no definido, ir a detalle para renombrar"
        return args

    def _file_posted(self):
        return len(request.files) > 0 or "filesAPI" in self.kwargs

    def _receive_and_prepare_file(self):
        import os
        formats = [".tif", ".gpkg", ".zip", ".json", ".geojson", ".csv", ".xlsx", ".xls"]
        path = ""
        folder = tempfile.mkdtemp(prefix=f"{app_acronym}_")
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
        try:
            gdf = gpd.read_file(path)
            return gdf
        except:
            try:
                gdf = import_pda_result(path, session=g.n_session.db_session)
                return gdf
            except:
                traceback.print_exc()
                return None

    def _read_raster_file(self, layer_name):
        from ..geo import geoserver_session
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


_ = register_api(bp_geo, LayersAPI, "geo/layers", f"{app_api_base}/geo/layers/", pk="_id", pk_type="string")
bp_geo.add_url_rule(app_api_base + '/geo/layers/<_id>.<string:_format>', view_func=_, methods=['GET'])


@Memoize
def get_styles():
    """ Dictionary of styles and a boolean telling if the palette is fixed or dynamic """
    palettes = {}
    # List of dynamic palettes
    try:
        sns.color_palette("", 3)
    except:
        s = traceback.format_exc()
        a = s.index('values are')
        s = s[a + 11:]
        a = s.index("\n")
        for palette in [_[1:] for _ in s[:a].split("', ")]:
            palettes[palette] = True
    # List of fixed palettes
    try:
        f_palettes = read_yaml(get_global_configuration_variable("COLOR_PALETTES"))
        for palette in f_palettes:
            palettes[palette["name"]] = False
    except:
        pass

    return palettes


def get_fixed_style(name):
    # palette_hex, intervals = get_fixed_style(color_ramp)
    # List of fixed palettes
    f_palettes = read_yaml(get_global_configuration_variable("COLOR_PALETTES"))
    colors = []
    intervals = []
    for palette in f_palettes:
        if name.lower() == palette["name"].lower():
            palette_type = palette.get("type", "continuous")
            for i in palette["intervals"]:
                colors.append(i["color"])
                if palette_type == "continuous":
                    intervals.append((i["left"], i["right"]))
                else:
                    intervals.append(str(i["value"]))  # CATEGORY read as String
    return colors, intervals


@n_session()
def get_styles_rest():
    # Return a list of available styles
    return ResponseObject(issues=[], status=200, content=get_styles(), count=0).get_response()


@n_session()
def put_property_attributes(id_, property_):
    from ..geo import geoserver_session
    issues = []
    status = 200
    req = request.get_json()
    palette = req.get("colormap")
    resp = dict()
    # Check if the layer "id_" and then, if the property exists
    session = g.n_session.db_session
    layer = session.query(GeographicLayer).filter(GeographicLayer.id == id_).first()
    property_dict = None
    if layer:
        for p_idx, p in enumerate(layer.properties):
            if p["name"] == property_:  # and p["data_type"] in ("numeric", "category"):
                property_dict = p
                break
    # Modify the property
    if property_dict:
        some_change = False
        for a in ["colormap", "data_type", "cat_type", "bins", "categories",
                  "col_type", "representable", "nature", "unit"]:
            if a in req:
                property_dict[a] = req[a]
                some_change = True
        if some_change:
            layer.properties[p_idx] = property_dict
            # Trick to mark "layer.properties" updatable for the session
            flag_modified(layer, "properties")

    # Check if the specified style exists
    if layer and property_dict and palette in get_styles():
        if p["data_type"] == "numeric" and not get_styles()[palette] or p["data_type"] == "category" and get_styles()[palette]:
            issues.append(Issue(IType.ERROR, f"The palette {palette} is not compatible with the property type {p['type']}"))
            status = 400
        else:
            style_name = f"layer_{id_}_{property_}"
            system_layer = True
            wkspc = workspace_names[0] if system_layer else workspace_names[1]
            geom_type = layer.attributes.get("geom_type", "Polygon") if layer.attributes else "Polygon"
            create_and_publish_style(geoserver_session,
                                     wkspc=wkspc,
                                     layer_name=layer.geoserver_name,
                                     attribute=property_,
                                     min_value=property_dict.get("min"),
                                     max_value=property_dict.get("max"),
                                     number_of_classes=req.get("bins", 7),
                                     categories=property_dict.get("categories"),
                                     style_name=style_name,
                                     color_ramp=palette,
                                     overwrite=True,
                                     geom_type=geom_type)
    else:
        if not layer:
            issues.append(Issue(IType.ERROR, f"Layer {id_} does not exist"))
        else:
            issues.append(Issue(IType.ERROR, f"Property {property_} does not exist or it is not 'numeric'"))
        if palette not in get_styles():
            issues.append(Issue(IType.ERROR, f"Palette {palette} does not exist"))

        status = 400

    return ResponseObject(issues=issues, status=status, content=resp, count=0).get_response()


bp_geo.add_url_rule(app_api_base + '/geo/layers/<int:id_>/<string:property_>',
                    view_func=put_property_attributes, methods=['PUT'])
bp_geo.add_url_rule(app_api_base + '/geo/layers/styles/',
                    endpoint="geo/layers/styles", view_func=get_styles_rest, methods=['GET'])


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
        from ..geo import geoserver_session
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
        from ..geo import geoserver_session
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
        from ..geo import geoserver_session
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
