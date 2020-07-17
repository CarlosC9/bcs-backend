from flask import Blueprint

bp_bos = Blueprint('api', __name__)


@bp_bos.route("/hello", methods=["GET"])
def hello():
    return "<h1 style='color:blue'>Hello hay!</h1>"


# from biobarcoding.services import SearchAPI
#
# search = SearchAPI.as_view('search_api')
#
# bp_bos.add_url_rule(
#     '/search',
#     view_func=search,
#     methods=['GET']
# )
# bp_bos.add_url_rule(
#     '/search/<string:type>',
#     view_func=search,
#     methods=['GET']
# )
