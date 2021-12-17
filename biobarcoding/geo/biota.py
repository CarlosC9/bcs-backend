import io
import json
import re
import traceback

import geopandas as gpd
import numpy as np
import pandas as pd

from ..services.species_names import get_canonical_species_names

"""
NOTE
GBIF service to parse species names, with 3 species from Biota:

curl "https://api.gbif.org/v1/parser/name?name=Euphorbia%20balsamifera%20Aiton%20subsp.%20balsamifera&name=Androcymbium%20hierrense%20A.%20Santos&name=Canaria%20tortuosa%20%28Webb%20%26%20Berthel.%29%20Jim.%20Mej%C3%ADas%20%26%20P.%20Vargas"
"""


def generate_pda_species_file_from_layer(db_sess, layer_id: int, layer_name: str, format_: str) -> str:
    from .. import postgis_engine, get_global_configuration_variable
    if not layer_name:
        return None

    if postgis_engine is None:
        from ..common.pg_helpers import create_pg_database_engine
        db_connection_string = get_global_configuration_variable('POSTGIS_CONNECTION_STRING')
        postgis_engine = create_pg_database_engine(db_connection_string, "ngd_geoserver", recreate_db=False)
    line_blocks = []
    if format_ == "nexus":
        line_blocks.append("#nexus\n\nbegin sets;")

    try:
        gdf = gpd.read_postgis(f"SELECT * FROM {layer_name} ORDER BY IDCELDA", postgis_engine, geom_col="geometry")
    except:
        return None

    if format_ == "species" or format_ == "species_canon":
        species = set()

    def adapt_encoding(s):
        try:
            return bytes(s.encode("iso-8859-1")).decode("utf8")
        except:
            return s

    # Map lower case column names to actual column names
    cols = {c.lower(): c for c in gdf.columns}
    # Check columns
    if "riqueza" in cols:
        underscores = format_ == "nexus"  # or format_ == "species_canon"
        # Unprocessed Biota. Read from "DENOMTAX" column
        for t, cell in gdf.groupby(cols["idcelda"]):
            if format_ == "nexus" or format_ == "pda_simple":
                # Re-encode
                _ = [adapt_encoding(taxon) for taxon in filter(None, cell[cols["denomtax"]].values)]
                # Normalize species names and discard None's (could not normalize species)
                _ = filter(None, get_canonical_species_names(db_sess, _, underscores))
                # Generate text lines for the cell ("celda")
                if format_ == "pda_simple":
                    line_blocks.append(str(len(cell)))
                    line_blocks.append('\n'.join(_))
                else:
                    id_celda = int(cell[cols['idcelda']].values[0])
                    _ = ' '.join(f"{sp}" for sp in _)
                    line_blocks.append(f"    taxset L{layer_id}_C{id_celda} = {_};")
            elif format_ == "species" or format_ == "species_canon":
                species.update([bytes(taxon.encode("iso-8859-1")).decode("utf8")
                                for taxon in filter(None, cell[cols["denomtax"]].values)
                                ])
    else:
        # Processed Biota. Read from "taxa" JSON column
        for t, r in gdf.iterrows():
            _ = json.loads(r["taxa"])
            if format_ == "nexus":
                pass
            elif format_ == "pda_simple":
                line_blocks.append(str(len(_)))
                line_blocks.append([r[cols["denomtax"]] for r in _.values()])

    if format_ == "nexus":
        line_blocks.append("end; [sets]")
    elif format_ == "species" or format_ == "species_canon":
        if format_ == "species":
            line_blocks.extend(species)
        else:
            line_blocks.extend(filter(None, get_canonical_species_names(db_sess, species)))
        line_blocks = sorted(line_blocks)

    return '\n'.join(line_blocks) + "\n", "text/plain"


# with open("/home/rnebot/Downloads/borrame.txt", "w") as f:
#     f.write(generate_pda_species_file_from_layer(layer_id=2, layer_name="layer_32", format_="pda_species"))


def import_pda_result(file_name: str, layer_id: str = None, session=None):
    # NOTE: this initial section is only for testing purposes
    from .. import engine, postgis_engine, get_global_configuration_variable
    from ..common.pg_helpers import create_pg_database_engine
    if engine is None:
        from ..db_models import DBSession
        db_connection_string = get_global_configuration_variable('DB_CONNECTION_STRING')
        engine = create_pg_database_engine(db_connection_string, "bcs", recreate_db=False)
        DBSession.configure(bind=engine)  # reconfigure the sessionmaker used by this scoped_session
        session = DBSession()
    if postgis_engine is None:
        db_connection_string = get_global_configuration_variable('POSTGIS_CONNECTION_STRING')
        postgis_engine = create_pg_database_engine(db_connection_string, "ngd_geoserver", recreate_db=False)

    # The function starts here

    # Read the file
    with open(file_name) as the_file:
        content = the_file.readlines()
    content = [x.strip() for x in content]
    equals_line = r"^\=+$"
    dashes_line: str = r"^\-+$"
    start_ = 0
    end_ = 0
    df = None
    while start_ < len(content):
        if re.match(equals_line, content[start_]):
            start_ += 2
            end_ = start_ + 1
            break
        start_ += 1
    if end_ < len(content):
        while end_ < len(content):
            if re.match(dashes_line, content[end_]):
                _ = io.StringIO()
                _.write('\n'.join(content[start_:end_]))
                _.seek(0)
                df = pd.read_csv(_, delim_whitespace=True)
                break
            end_ += 1
    if df is not None:
        df.columns = ["name", "pd", "pd_excl", "pd_endem"]
        # Find layer ID
        if layer_id is None:
            tmp = re.match(r"L([0-9]+)_C([0-9]+)", df.iloc[0, 0])
            if tmp:
                layer_id = int(tmp.groups()[0])
        if layer_id is None:
            return None  # Layer ID not present as parameter or in the input file
        else:
            # Find layer_name
            from ..db_models.geographics import GeographicLayer
            layer = session.query(GeographicLayer).filter(GeographicLayer.id == layer_id).first()
            layer_name = layer.geoserver_name
            # Read reference Biota layer, to merge geometry
            try:
                gdf = gpd.read_postgis(f"SELECT DISTINCT(idcelda), geometry FROM {layer_name} ORDER BY IDCELDA",
                                       postgis_engine, geom_col="geometry")
                gdf["idcelda"] = gdf["idcelda"].astype(np.int64)
            except:
                traceback.print_exc()
                return None  # Could not read the source Biota layer

            # Generate column "idcelda" from "name" column
            df["idcelda"] = df["name"].str[len(f"L{layer_id}_") + 1:].astype(np.int64)
            del df["name"]
            # Join and return the layer
            return gdf.merge(df, on="idcelda")
    else:
        return None  # Could not parse the PDA output file


# import_pda_result("/home/rnebot/GoogleDrive/AA_NEXTGENDEM/pda.out")


# TODO Specific to BIOTA, should not be in the generic GeoAPI but in a special preprocessing module
def read_biota_file(path):
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
