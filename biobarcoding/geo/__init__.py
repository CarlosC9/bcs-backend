from biobarcoding import postgis_engine
from sqlalchemy import create_engine
from geo.Geoserver import Geoserver
import geopandas as gpd

'''
geo_sesion inizialize (just object with admin and password information)

'''
geoserver_session = None

layers = {"plantas":"/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/biobarcoding/geo_layers/plantas/Plantas.shp"}

def inizialice_layers(app):
    engine = create_engine(app.config['POSTGIS_CONNECTION_STRING'] + "ngd_geoserver")
    for layer_name, path in layers.items():
        df = gpd.read_file(path)
        try:
            df.to_postgis(layer_name,engine, if_exists="fail")
        except ValueError as v:
            print(v)
            pass


def inizialice_geoserver(flask_app):
    global geoserver_session
    inizialice_layers(flask_app)
    if {'GEOSERVER_USER',
        'GEOSERVER_PASSWORD',
        'GEOSERVER_HOST',
        'GEOSERVER_PORT'} <= flask_app.config.keys():
        geoserver_url = f"http://{flask_app.config['GEOSERVER_HOST']}:{flask_app.config['GEOSERVER_PORT']}/geoserver"
        geoserver_session = Geoserver(geoserver_url, username= flask_app.config['GEOSERVER_USER'], password= flask_app.config['GEOSERVER_PASSWORD'])
        workspaces = geoserver_session.get_workspaces()
        print(workspaces)
        if workspaces.get('workspaces') != '': #there at list a workspace
            if {'name': 'ngd', 'href': f'{geoserver_url}/rest/workspaces/ngd.json'} in workspaces['workspaces'].get('workspace'):
                print('workspace ngd ready')
        else:
            geoserver_session.create_workspace(workspace='ngd')
            # the storage is always asociated to the workspace so, no need to ask if the feature storage already exists
            # TODO ask for any other storage existance?
            geoserver_session.create_featurestore(store_name = 'geo_data',
                                    workspace='ngd',
                                    db = flask_app.config['POSTGIS_DB'],
                                    host = 'postgis', # internal host TODO get this as variable
                                    port = 5435, # internal port TODO get this as variable
                                    pg_user= flask_app.config['POSTGIS_USER'],
                                    pg_password= flask_app.config['POSTGIS_PASSWORD']
                                )
        # todo check that all initial layers are publish
        # layers_in_geoserver = geo.get_layers()

        for layer_name, _ in layers.items():
            # todo error al publicar
            geoserver_session.publish_featurestore(workspace='ngd', store_name='geo_data', pg_table=layer_name)
    else:
        print("no geoserver data in config file cant open GEOSERVER session")

