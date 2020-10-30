from flask import Blueprint, abort, request, make_response, jsonify, session as flask_session
from flask.views import MethodView

from biobarcoding.authentication import obtain_identity_from_request, serialize_session, BCSSession, deserialize_session
from biobarcoding.rest import bcs_api_base, register_api
from firebase_admin import auth

bp_auth = Blueprint('bp_auth', __name__)


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
    response_object = {
        'status': 'success',
        'message': 'valid token'
    }
    return make_response(jsonify(response_object)), 200


class AuthnAPI(MethodView):
    """
    Authentication service. In charge of creating initial BCSSession, so "bcs_session" decorator can
    use it to check authorization of functions
    """

    decorators = []  # Add decorators: identity, function execution permissions, logging, etc.

    def get(self):
        """ Obtain current Identity or 'anonymous' """
        sess = deserialize_session(flask_session.get("session"), False)
        if sess is not None and isinstance(sess, BCSSession):
            # TODO Return current identity
            response_object = {
                'status': 'success',
                'message': f"Identity: {sess.identity_name}"
            }
        else:
            # TODO Return "not logged in"
            response_object = {
                'status': 'success',
                'message': ""
            }
        return make_response(jsonify(response_object)), 200

    def post(self):
        """ Which sense does this have? """
        response_object = {
            'status': 'error',
            'message': 'not implemented'
        }

        return make_response(jsonify(response_object)), 200

    def delete(self):
        """ Logout """
        flask_session["session"] = None
        response_object = {
            'status': 'success',
            'message': ""
        }
        return make_response(jsonify(response_object)), 200

    def put(self):
        """ Login """
        identity = obtain_identity_from_request()
        # TODO If identity does not exist, create one, and a relation with the authentication method
        sess = BCSSession()
        sess.identity_name = identity
        flask_session["session"] = serialize_session(sess)
        # Attach identity to the current session
        response_object = {
            'status': 'success',
            'message': ""
        }
        return make_response(jsonify(response_object)), 200


# Special behavior: "authn" is a singleton, which can be None or defined with a login
# POST not implemented
view_func = AuthnAPI.as_view("authn")
bp_auth.add_url_rule(f"{bcs_api_base}/authn", view_func=view_func, methods=['GET', 'PUT', 'DELETE'])

