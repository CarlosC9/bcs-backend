from flask import (Flask, request, session as flask_session, redirect, current_app)
from flask_session import Session as FlaskSessionServerSide
from flask_cors import CORS

import biobarcoding
from biobarcoding.jobs.galaxy_resource import initialize_galaxy
from biobarcoding.rest import logger, log_level, load_configuration_file, construct_session_persistence_backend, \
    initialize_database, initialize_database_chado, bcs_gui_base
from biobarcoding.rest.auth import bp_auth
from biobarcoding.rest.identities_and_company import bp_identities, bp_sys_functions, bp_roles, bp_identities_roles, \
    bp_groups, bp_organizations, bp_acl
from biobarcoding.rest.sequences import bp_sequences
from biobarcoding.rest.alignments import bp_alignments
from biobarcoding.rest.phylotrees import bp_phylotrees
from biobarcoding.rest.taxonomies import bp_taxonomies
from biobarcoding.rest.organisms import bp_organisms
from biobarcoding.rest.analyses import bp_analyses
from biobarcoding.rest.ontologies import bp_ontologies
from biobarcoding.rest.jobs import bp_jobs
from biobarcoding.rest.tasks import bp_tasks
from biobarcoding.rest.processes import  bp_processes, bp_resources
from biobarcoding.rest.gui_static import bp_gui
from biobarcoding.tasks import initialize_celery
from biobarcoding.authentication import initialize_firebase

# Flask and configuration file
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
    CORS(app,
         resources={r"/api/*": {"origins": "*"}},
         supports_credentials=True
         )

    # Database BCS
    initialize_database(app)

    # Database Chado
    initialize_database_chado(app)

    # Security
    # initialize_authn_authr(app)

    # RESTful endpoints
    for bp in [bp_auth,
               bp_sequences,
               bp_alignments,
               bp_phylotrees,
               bp_taxonomies,
               bp_organisms,
               bp_analyses,
               bp_ontologies,
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
               bp_resources
               ]:
        app.register_blueprint(bp)

    # Celery
    initialize_celery(app)

    # Galaxy
    print("Initializing base Galaxy instance")
    initialize_galaxy(app)
    print("Initializing base Galaxy instance - DONE")

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
    from biobarcoding.tasks.definitions import add
    # add.delay(2, 3)
    return f"<h1 style='color:blue'>{biobarcoding.flask_app.config['DB_CONNECTION_STRING']}!<br>{current_app.config['DB_CONNECTION_STRING']}</h1>"


@biobarcoding.flask_app.after_request
def after_a_request(response):
    # Keep cookies ...
    for i in request.cookies.items():
        response.set_cookie(i[0], i[1])

    # ... except when Invalidation requested
    if "__invalidate__" in flask_session:
        response.delete_cookie(current_app.session_cookie_name)

    return response


if __name__ == "__main__":
    biobarcoding.flask_app.run(host='0.0.0.0')
