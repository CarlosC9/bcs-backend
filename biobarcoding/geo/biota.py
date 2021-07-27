import json
import geopandas as gpd
import numpy as np


def generate_pda_species_file_from_layer(layer_name: str) -> str:
    from biobarcoding import postgis_engine, get_global_configuration_variable
    if not layer_name:
        return None

    if postgis_engine is None:
        from biobarcoding.common.pg_helpers import create_pg_database_engine
        db_connection_string = get_global_configuration_variable('POSTGIS_CONNECTION_STRING')
        postgis_engine = create_pg_database_engine(db_connection_string, "ngd_geoserver", recreate_db=False)
    line_blocks = []
    gdf = gpd.read_postgis(f"SELECT * FROM {layer_name} ORDER BY IDCELDA", postgis_engine, geom_col="geometry")
    # Map lower case column names to actual column names
    cols = {c.lower(): c for c in gdf.columns}
    # Check columns
    if "riqueza" in cols:
        # Unprocessed Biota. Read from "DENOMTAX" column
        for t, cell in gdf.groupby(cols["idcelda"]):
            line_blocks.append(str(len(cell)))
            line_blocks.append('\n'.join(filter(None, cell[cols["denomtax"]].values)))
    else:
        # Processed Biota. Read from "taxa" JSON column
        for t, r in gdf.iterrows():
            _ = json.loads(r["taxa"])
            line_blocks.append(str(len(_)))
            line_blocks.append([r[cols["denomtax"]] for r in _.values()])

    return '\n'.join(line_blocks)


# with open("/home/rnebot/Downloads/borrame.txt", "w") as f:
#     f.write(generate_pda_species_file_from_layer(layer_name="layer_32"))


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
