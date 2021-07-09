from geo.Geoserver import Geoserver


'''
https://geoserver-rest.readthedocs.io/en/latest/how_to_use.html#getting-started-with-geoserver-rest
geo_sesion inizialization (just object with admin and password information)

IMPORTANT!:
Expose Primary Key in geoserserver 
http://localhost:9180/geoserver/web/wicket/bookmarkable/org.geoserver.web.data.store.DataAccessEditPage?8&storeName=geo_data&wsName=ngd
check Expose primary keys
'''
geoserver_session = None


def initialize_geoserver(flask_app):
    global geoserver_session
    if {'GEOSERVER_USER',
        'GEOSERVER_PASSWORD',
        'GEOSERVER_HOST',
        'GEOSERVER_PORT'} <= flask_app.config.keys():
        geoserver_url = f"http://{flask_app.config['GEOSERVER_HOST']}:{flask_app.config['GEOSERVER_PORT']}/geoserver"
        geoserver_session = Geoserver(geoserver_url, username=flask_app.config['GEOSERVER_USER'], password=flask_app.config['GEOSERVER_PASSWORD'])
        workspaces = geoserver_session.get_workspaces()
        print(workspaces)
        # CHECK GEOSERVER:
        # get geoserver version
        version = geoserver_session.get_version()
        print(version)
        # # get ststem info
        # status = geoserver_session.get_status()
        # print(status)
        # system_status = geoserver_session.get_system_status()
        # print(system_status)
        if workspaces.get('workspaces') != '' and \
                {'name': 'ngd', 'href': f'{geoserver_url}/rest/workspaces/ngd.json'} in workspaces['workspaces'].get('workspace'):
            print('Workspace "ngd" ready')
        else:
            geoserver_session.create_workspace(workspace='ngd')

        ds = geoserver_session.get_datastores('ngd')
        if ds.get('dataStores') == '':
            # the storage is always associated to the workspace so, no need to ask if the feature storage already exists
            geoserver_session.create_featurestore(store_name='geo_data',
                                                  workspace='ngd',
                                                  db=flask_app.config['POSTGIS_DB'],
                                                  host=flask_app.config['POSTGIS_HOST'],
                                                  port=int(flask_app.config['POSTGIS_PORT']),
                                                  pg_user=flask_app.config['POSTGIS_USER'],
                                                  pg_password=flask_app.config['POSTGIS_PASSWORD']
                                                  )
    else:
        print("no geoserver data in config file cant open GEOSERVER session")

