from flask import (Flask, Response, request, session as flask_session, send_from_directory
                   )
from flask_session import Session as FlaskSessionServerSide
from flask_cors import CORS

from biobarcoding.rest import logger, log_level, load_configuration_file, construct_session_persistence_backend, \
    initialize_database
from biobarcoding.rest.bos import bp_bos
from biobarcoding.rest.gui_static import bp_gui


app = Flask(__name__)
app.debug = True
UPLOAD_FOLDER = '/tmp/'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
load_configuration_file(app)
initialize_database(app)
d = construct_session_persistence_backend(app)
app.config.update(d)

FlaskSessionServerSide(app)  # Flask Session
CORS(app,                    # CORS
     resources={r"/api/*": {"origins": "*"}},
     supports_credentials=True
     )
app.logger.setLevel(log_level)
logger.setLevel(log_level)

for bp in [bp_bos, bp_gui]:
    app.register_blueprint(bp)


@app.route("/test")
def test_rest_open():
    return "<h1 style='color:blue'>Test!</h1>"


@app.after_request
def after_a_request(response):
    # Keep cookies ...
    for i in request.cookies.items():
        response.set_cookie(i[0], i[1])

    # ... except when Invalidation requested
    if "__invalidate__" in flask_session:
        response.delete_cookie(app.session_cookie_name)

    return response


if __name__ == "__main__":
    app.run(host='0.0.0.0')
