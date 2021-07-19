from geo.Geoserver import Geoserver


'''
https://geoserver-rest.readthedocs.io/en/latest/how_to_use.html#getting-started-with-geoserver-rest
geo_sesion inizialization (just object with admin and password information)

IMPORTANT!:
Expose Primary Key in geoserver 
http://localhost:9180/geoserver/web/wicket/bookmarkable/org.geoserver.web.data.store.DataAccessEditPage?8&storeName=geo_data&wsName=ngd
check Expose primary keys
'''
geoserver_session = None
workspace_names = ["ngd", "ngd_users"]
postgis_store_name = f"postgis"


def initialize_geoserver(flask_app):
    global geoserver_session
    if {'GEOSERVER_USER',
        'GEOSERVER_PASSWORD',
        'GEOSERVER_HOST',
        'GEOSERVER_PORT'} <= flask_app.config.keys():
        geoserver_url = f"http://{flask_app.config['GEOSERVER_HOST']}:{flask_app.config['GEOSERVER_PORT']}/geoserver"
        geoserver_session = Geoserver(geoserver_url, username=flask_app.config['GEOSERVER_USER'], password=flask_app.config['GEOSERVER_PASSWORD'])
        # CHECK GEOSERVER:
        # get geoserver version
        version = geoserver_session.get_version()
        print(version)
        # Get system info
        status = geoserver_session.get_status()
        print(status)
        system_status = geoserver_session.get_system_status()
        print(system_status)
        # Create workspaces
        workspaces = geoserver_session.get_workspaces()
        for wkspc in workspace_names:
            if workspaces.get('workspaces') != '' and \
                    {'name': wkspc, 'href': f'{geoserver_url}/rest/workspaces/{wkspc}.json'} in workspaces['workspaces'].get('workspace'):
                print(f'Workspace {wkspc} ready')
            else:
                geoserver_session.create_workspace(workspace=wkspc)

            ds = geoserver_session.get_datastores(wkspc)
            if ds.get('dataStores') == '':
                # The storage is always associated to the workspace so,
                # no need to ask if the feature storage already exists
                geoserver_session.create_featurestore(store_name=postgis_store_name,
                                                      workspace=wkspc,
                                                      db=flask_app.config['POSTGIS_DB'],
                                                      host=flask_app.config['POSTGIS_HOST'],
                                                      port=int(flask_app.config['POSTGIS_PORT']),
                                                      pg_user=flask_app.config['POSTGIS_USER'],
                                                      pg_password=flask_app.config['POSTGIS_PASSWORD']
                                                      )
    else:
        print("no geoserver data in config file cant open GEOSERVER session")

