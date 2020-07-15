from flask import Blueprint

bp_bos = Blueprint('api', __name__)


@bp_bos.route("/hello", methods=["GET"])
def hello():
    return "<h1 style='color:blue'>Hello hay!</h1>"
