from pathlib import Path
import sys

from NamedAtomicLock import NamedAtomicLock
from flask_session import Session as FlaskSessionServerSide
from flask import (Flask, request, session as flask_session, redirect, current_app)
from flask_cors import CORS
from flask_socketio import SocketIO

# Workaround for relative imports in a "__main__"
# From: https://stackoverflow.com/a/28154841

if __name__ == '__main__' and (__package__ is None or __package__ == ""):
    file = Path(__file__).resolve()
    parent, top = file.parent, file.parents[3]

    sys.path.append(str(top))
    try:
        sys.path.remove(str(parent))
    except ValueError:  # Already removed
        pass

    import biobarcoding.rest
    __package__ = 'biobarcoding.rest'
# End of workaround

import biobarcoding as base_app_pkg
from ..authentication import initialize_firebase
from ..geo import initialize_geoserver
from ..rest.views_dashboards import bp_viewz, bp_dashboards, bp_dashboard_items
from . import logger, log_level, load_configuration_file, construct_session_persistence_backend, initialize_database, app_gui_base, ResponseObject, init_socket, initialize_postgis, initialize_ssh, \
    initialize_chado_edam, initialize_database_chado, initialize_galaxy
from biobarcoding import app_acronym
from ..services.geoprocesses import update_geoprocesses
from .auth import bp_auth, bp_api_key
from .browser_filters import bp_bfilters
from .files import bp_files
from .geo_rest import bp_geo
from .geoprocesses_and_case_studies import bp_case_studies, bp_geoprocesses, bp_geoprocesses_ports, bp_geoprocess_instances, \
    bp_geoprocess_port_types, bp_case_studies_fos
from .gui_static import bp_gui
from .hierarchies import bp_hierarchies, bp_hierarchy_nodes
from .identities_and_company import bp_identities, bp_sys_functions, bp_roles, bp_identities_roles, \
    bp_groups, bp_organizations, bp_acl, bp_identities_authenticators, bp_identity_store
from .jobs import bp_jobs
from .bio import bp_bio
from .annotation_forms import bp_annotations
from .processes import bp_processes, bp_resources
from .proxy import bp_proxy
from .tasks import bp_tasks
from .views import bp_views
from ..tasks import initialize_celery


# Flask and configuration file
socket_service_socketio = None
app = None


def create_app(debug, start_socket=True, cfg_dict=None):
    """

    TODO - Possible improvements, noted at: http://www.patricksoftwareblog.com/structuring-a-flask-project/

    :return:
    """

    global app
    app = Flask(__name__)
    app.debug = debug

    UPLOAD_FOLDER = '/tmp/'
    app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
    load_configuration_file(app)

    initialize_firebase(app)

    # Session persistence configuration
    d = construct_session_persistence_backend(app)
    app.config.update(d)

    # Flask Session
    FlaskSessionServerSide(app)

    # CORS
    CORS(app,
         resources={r"/api/*": {"origins": "*"}, r"/pxy/*": {"origins": "*"}},
         supports_credentials=True
         )

    if start_socket:
        global socket_service_socketio
        socket_service_socketio = SocketIO(app, cors_allowed_origins="*", async_mode="threading")
        init_socket(socket_service_socketio)

    lock = NamedAtomicLock(f"{app_acronym}-backend-lock")
    lock.acquire()
    try:
        # Database
        initialize_database(app)

        # Database Chado
        initialize_database_chado(app)

        # Insert EDAM Ontology in chado
        initialize_chado_edam(app)

        # Insert BibTeX annotation forms
        from biobarcoding.forms import initialize_bibtex_forms
        initialize_bibtex_forms()

        # PostGIS Database
        initialize_postgis(app)

        # Galaxy
        print("Initializing Galaxy instances")
        # initialize_galaxy(app)
        print("Initializing Galaxy instances - DONE")

        # SSH
        print("Initializing SSH resources connections")
        initialize_ssh(app)
        print("Initializing SSH resources connections - DONE")
    finally:
        lock.release()

    # Following code was used (modify layer ID) to research ISSUE of extra slow export in some cases
    # The reason was a version of PyProj (3.1.0 slow, 2.6.1 normal -fast-)
    #
    # from .geo_rest import LayersAPI, get_content
    # from base_app_pkg.db_models import DBSession
    # from base_app_pkg.db_models.geographics import GeographicLayer
    # issues = []
    # issues, layer, count, status = get_content(DBSession(), GeographicLayer, issues, "14")
    # content = LayersAPI()._export(DBSession(), layer, _format="nexus")
    # sys.exit(1)

    print("Initializing Geoprocesses - requesting info to geoprocess executors")
    update_geoprocesses()
    print("Initializing Geoprocesses - DONE")

    print("Initializing Geoserver")
    initialize_geoserver(app)
    print("Initializing Geoserver - DONE")

    # Security
    # initialize_authn_authr(app)

    # RESTful endpoints
    for bp in [
               bp_auth,
               bp_files,
               bp_jobs,
               bp_gui,
               bp_hierarchies,
               bp_hierarchy_nodes,
               bp_identities,
               bp_identities_authenticators,
               bp_sys_functions,
               bp_roles,
               bp_identities_roles,
               bp_groups,
               bp_organizations,
               bp_acl,
               bp_processes,
               bp_resources,
               bp_annotations,
               bp_bfilters,
               bp_geo,
               bp_views,
               bp_proxy,
               bp_case_studies,
               bp_case_studies_fos,
               bp_geoprocesses,
               bp_geoprocesses_ports,
               bp_geoprocess_instances,
               bp_geoprocess_port_types,
			   bp_bio,
			   bp_tasks,
               bp_identity_store,
               bp_viewz,
               bp_dashboards,
               bp_dashboard_items,
               bp_api_key,
			   ]:
        app.register_blueprint(bp)

    # Celery
    initialize_celery(app)

    # Logger
    app.logger.setLevel(log_level)

    return app


# FLASK_ENV=development FLASK_APP=biobarcoding.rest.main flask run
logger.debug("BACKEND is STARTING...")

start_socket = True
base_app_pkg.flask_app = create_app(True, start_socket=start_socket)

logger.debug("BACKEND INITIALIZATION FINISHED !!!")


@base_app_pkg.flask_app.route("/")
def index():
    return redirect(app_gui_base)


@base_app_pkg.flask_app.route("/test")
def test_rest_open():
    r = ResponseObject()
    r.content = dict(it_works="yes!")
    return r.get_response()


@base_app_pkg.flask_app.after_request
def after_a_request(response):
    return response

    # Keep cookies ...
    for i in request.cookies.items():
        response.set_cookie(i[0], i[1])

    # ... except when Invalidation requested
    if "__invalidate__" in flask_session:
        response.delete_cookie(current_app.session_cookie_name)
    else:
        # Allow Cross Site usage when debugging
        if base_app_pkg.get_global_configuration_variable("SAMESITE_NONE", "False") == "True":
            for i, h in enumerate(response.headers):
                if h[0].lower() == "set-cookie" and h[1].startswith(f"{current_app.session_cookie_name}="):
                    response.headers[i] = (h[0], f"{h[1]}; SameSite=None")

    return response


if __name__ == "__main__":
    localhost_name = "0.0.0.0"  # "localhost"  # 0.0.0.0, 127.0.0.1
    if not start_socket:
        logger.debug("BACKEND is STARTING FLASK SERVER !!!")
        base_app_pkg.flask_app.run(host=localhost_name,
                                   use_reloader=False,  # Avoid loading twice the application
                                   )
    else:
        logger.debug("BACKEND is STARTING FLASK SERVER WITH WEBSOCKETS !!!")
        socket_service_socketio.run(base_app_pkg.flask_app, host=localhost_name, use_reloader=False)
