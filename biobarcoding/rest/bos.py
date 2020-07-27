from flask import Blueprint

bp_bos = Blueprint('api', __name__)


@bp_bos.route("/hello", methods=["GET"])
def hello():
    return "<h1 style='color:blue'>Hello hay!</h1>"


# from flask import request, make_response, jsonify
# from flask.views import MethodView
#
# class SearchAPI(MethodView):
#     """
#     Search Resource
#     """
#     def get(self, type=None):
#         msg = f'GET {request.path}\nGetting search {type}'
#         print(msg)
#         self._check_data()
#
#         responseObject = {
#             'status': 'success',
#             'message': msg
#         }
#         return make_response(jsonify(responseObject)), 200
#
#
#     def _check_data(self):
#         post_data = request.get_json()
#         print(f'JSON data: {post_data}')
#

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
