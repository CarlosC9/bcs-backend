from flask import (Flask, request, session as flask_session, redirect, current_app)
from flask_session import Session as FlaskSessionServerSide
from flask_cors import CORS
from NamedAtomicLock import NamedAtomicLock
from flask_socketio import SocketIO

import biobarcoding
from biobarcoding.rest import logger, log_level, load_configuration_file, construct_session_persistence_backend, \
    initialize_database, initialize_database_chado, bcs_gui_base, ResponseObject, initialize_chado_edam, \
    init_socket, initialize_postgis, initialize_ssh, initialize_galaxy
from biobarcoding.rest.auth import bp_auth
from biobarcoding.rest.file_manager import bp_files
from biobarcoding.rest.identities_and_company import bp_identities, bp_sys_functions, bp_roles, bp_identities_roles, \
    bp_groups, bp_organizations, bp_acl
from biobarcoding.rest.bos import bp_bos
from biobarcoding.rest.metadata import bp_metadata
from biobarcoding.rest.jobs import bp_jobs
from biobarcoding.rest.tasks import bp_tasks
from biobarcoding.rest.processes import  bp_processes, bp_resources
from biobarcoding.rest.geo_rest import bp_geo
from biobarcoding.rest.gui_static import bp_gui
from biobarcoding.rest.browser_filters import bp_bfilters
from biobarcoding.rest.views import bp_views
from biobarcoding.rest.proxy import bp_proxy
from biobarcoding.tasks import initialize_celery
from biobarcoding.authentication import initialize_firebase
from biobarcoding.geo import initialize_geoserver

# Flask and configuration file
socket_service_socketio = None
app = None


def create_app(debug, cfg_dict=None):
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
    cors = CORS(app,
         resources={r"/api/*": {"origins": "*"}},
         supports_credentials=True
         )

    global socket_service_socketio
    socket_service_socketio = SocketIO(app, cors_allowed_origins="*", async_mode="threading")
    init_socket(socket_service_socketio)

    lock = NamedAtomicLock("bcs-backend-lock")
    lock.acquire()
    try:
        # Database BCS
        initialize_database(app)

        # Database Chado
        initialize_database_chado(app)

        # Insert EDAM Ontology in chado
        initialize_chado_edam(app)

        # iniitalize postgis Database
        initialize_postgis(app)

        # Galaxy
        print("Initializing Galaxy instances")
        initialize_galaxy(app)
        print("Initializing Galaxy instances - DONE")

        # SSH
        print("Initializing SSH resources connections")
        initialize_ssh(app)
        print("Initializing SSH resources connections - DONE")

    finally:
        lock.release()

    print("Initializing Geoserver")
    initialize_geoserver(app)
    print("Initializing Geoserver - DONE")

    # Security
    # initialize_authn_authr(app)

    # RESTful endpoints
    for bp in [bp_auth,
               bp_files,
               bp_jobs,
               bp_tasks,
               bp_gui,
               bp_identities,
               bp_sys_functions,
               bp_roles,
               bp_identities_roles,
               bp_groups,
               bp_organizations,
               bp_acl,
               bp_processes,
               bp_resources,
               bp_bos,
               bp_metadata,
               bp_bfilters,
               bp_geo,
               bp_views,
               bp_proxy,
               ]:
        app.register_blueprint(bp)

    # Celery
    initialize_celery(app)

    # Logger
    app.logger.setLevel(log_level)
    logger.setLevel(log_level)

    return app


biobarcoding.flask_app = create_app(True)


@biobarcoding.flask_app.route("/")
def index():
    return redirect(bcs_gui_base)


@biobarcoding.flask_app.route("/test")
def test_rest_open():
    r = ResponseObject()
    r.content = dict(it_works="yes!", other_info=dict(conn=biobarcoding.flask_app.config['DB_CONNECTION_STRING']))
    return r.get_response()


@biobarcoding.flask_app.after_request
def after_a_request(response):
    # Keep cookies ...
    for i in request.cookies.items():
        response.set_cookie(i[0], i[1])

    # ... except when Invalidation requested
    if "__invalidate__" in flask_session:
        response.delete_cookie(current_app.session_cookie_name)
    else:
        # Allow Cross Site usage when debugging
        if biobarcoding.get_global_configuration_variable("SAMESITE_NONE", "False") == "True":
            for i, h in enumerate(response.headers):
                if h[0].lower() == "set-cookie" and h[1].startswith(f"{current_app.session_cookie_name}="):
                    response.headers[i] = (h[0], f"{h[1]}; SameSite=None")

    return response

if __name__ == "__main__":
    biobarcoding.flask_app.run(host='0.0.0.0',
                               use_reloader=False,  # Avoid loading twice the application
                               )
    socket_service_socketio.run(biobarcoding.flask_app, host='0.0.0.0')
