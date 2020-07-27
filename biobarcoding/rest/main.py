from flask import (Flask, Response, request, session as flask_session, send_from_directory
                   )
from flask_session import Session as FlaskSessionServerSide
from flask_cors import CORS

from biobarcoding.rest import logger, log_level, load_configuration_file, construct_session_persistence_backend, \
    initialize_database
from biobarcoding.rest.auth import bp_auth
from biobarcoding.rest.bos import bp_bos
from biobarcoding.rest.seq import bp_seq
from biobarcoding.rest.msa import bp_msa
from biobarcoding.rest.phylo import bp_phylo
from biobarcoding.rest.ansis import bp_ansis
from biobarcoding.rest.onto import bp_onto
from biobarcoding.rest.taxon import bp_taxon
from biobarcoding.rest.io import bp_io
from biobarcoding.rest.job import bp_job
from biobarcoding.rest.jobqueue import bp_jobqueue
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

# for bp in [bp_bos, bp_gui]:
for bp in [bp_auth,bp_bos,bp_seq,bp_msa,bp_phylo,bp_ansis,bp_onto,bp_taxon,bp_io,bp_job,bp_jobqueue,bp_gui]:
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
