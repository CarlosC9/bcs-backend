from flask import (Flask, request, session as flask_session, redirect, current_app)
from flask_session import Session as FlaskSessionServerSide
from flask_cors import CORS

import biobarcoding
from biobarcoding.rest import logger, log_level, load_configuration_file, construct_session_persistence_backend, \
    initialize_database, bcs_gui_base
from biobarcoding.rest.auth import bp_auth
from biobarcoding.rest.bos import bp_bos
from biobarcoding.rest.seq import bp_seq
from biobarcoding.rest.msa import bp_msa
from biobarcoding.rest.phylo import bp_phylo
from biobarcoding.rest.meta import bp_meta
from biobarcoding.rest.in_out import bp_io
from biobarcoding.rest.jobs import bp_jobs
from biobarcoding.rest.tasks import bp_tasks
from biobarcoding.rest.gui_static import bp_gui
from biobarcoding.tasks import initialize_celery
from biobarcoding.authentication import initialize_firebase

# Flask and configuration file


def create_app(debug, cfg_dict=None):
    """

    TODO - Possible improvements, noted at: http://www.patricksoftwareblog.com/structuring-a-flask-project/

    :return:
    """
    app = Flask(__name__)
    app.debug = debug

    UPLOAD_FOLDER = '/tmp/'
    app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

    load_configuration_file(app)

    initialize_firebase(app)

    FlaskSessionServerSide(app)  # Flask Session
    CORS(app,                    # CORS
         resources={r"/api/*": {"origins": "*"}},
         supports_credentials=True
         )

    for bp in [bp_auth,
               bp_bos,
               bp_seq,
               bp_msa,
               bp_phylo,
               bp_meta,
               bp_io,
               bp_jobs,
               bp_tasks,
               bp_gui
               ]:
        app.register_blueprint(bp)

    # Database
    initialize_database(app)

    # Session persistence
    d = construct_session_persistence_backend(app)
    app.config.update(d)

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
