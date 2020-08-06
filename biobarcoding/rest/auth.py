from flask import Blueprint, abort, request, make_response, jsonify
from biobarcoding.rest import bcs_api_base

bp_auth = Blueprint('auth', __name__)

from firebase_admin import auth

# Access from bcs-sys through 'rev-proxy/conf.d/sub-auth.conf'
@bp_auth.route("/auth", methods=["GET","POST","PUT","DELETE"])
def token_verification():
    auth_token = None
    if 'Authorization' in request.headers:
        auth_token = request.headers['Authorization'].split(" ")[1]
    elif 'X-API-Key' in request.headers:
        abort(405, 'The X-API-Key authentication is not available.')
    if not auth_token:
        abort(401, 'The session token is missing')
    try:
        print(auth.verify_id_token(auth_token))
    except Exception as e:
        print(e)
        abort(401, 'The session token is not valid or has expired')
    responseObject = {
        'status': 'success',
        'message': 'valid token'
    }
    return make_response(jsonify(responseObject)), 200


# from flask.views import MethodView
#
# class AuthAPI(MethodView):
#     """
#     Authentication Resource
#     """
#     def post(self):
#         msg = f'POST {request.path}\nLogging in or out.'
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
#         if request.path == '/auth/login':
#             msg = 'Successfully logged in.'
#         elif request.path == '/auth/logout':
#             msg = 'Successfully logged in.'
#         else:
#             msg = 'Warning: Route may be wrong.'
#         try:
#             post_data = request.get_json()
#             print(f'{msg}\nJSON data: {post_data}')
#         except Exception as e:
#             pass
#
#
# auth = AuthAPI.as_view('auth_api')
# bp_auth.add_url_rule(
#     bcs_api_base + '/auth/login',
#     view_func=auth,
#     methods=['POST']
# )
# bp_auth.add_url_rule(
#     bcs_api_base + '/auth/logout',
#     view_func=auth,
#     methods=['POST']
# )
