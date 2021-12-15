from datetime import datetime

from firebase_admin import auth
from flask import Blueprint, abort, request, make_response, jsonify, session as flask_session, g
from flask.views import MethodView
from sqlalchemy import and_, or_

from ..authentication import serialize_session, AppSession, deserialize_session, obtain_idauth_from_request, n_session
from . import app_api_base, ResponseObject, logger
from ..authorization import string_to_ast, authr_expression, ast_evaluator
from ..common.helpers import obj_to_json
from ..db_models import ObjectType
from ..db_models.sysadmin import SystemFunction, ACLExpression, ACL

bp_auth = Blueprint('bp_auth', __name__)


# Access from <app>-sys through '<app>_ecosys/conf.d/sub-auth.conf'
@bp_auth.route("/auth", methods=["GET", "POST", "PUT", "DELETE"])
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
    Authentication service. In charge of creating initial AppSession, so "n_session" decorator can
    use it to check authorization of functions
    """

    decorators = []  # Add decorators: identity, function execution permissions, logging, etc.

    def get(self):
        """ Obtain current Identity or 'anonymous' """
        sess = deserialize_session(flask_session.get("session"), False)
        if sess is not None and isinstance(sess, AppSession):
            # TODO Return current identity
            response_object = {
                'status': 'success',
                'message': f"Identity: {sess.identity_name}",
                'identity': sess.identity_name,
                'identity_id': sess.identity_id
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
        # If identity does not exist, create one, and a relation with the authentication method
        id_auth = obtain_idauth_from_request()  # MAIN !!
        # If the identity is Active, continue;
        # If not, return an error
        if id_auth:
            # logger.debug(f"AUTH: {obj_to_json(id_auth)}")
            # logger.debug(f"IDENTITY: {obj_to_json(id_auth.identity)}")
            if not id_auth.identity.deactivation_time:  # and id_auth.identity.can_login
                sess = AppSession()
                sess.identity_id = id_auth.identity.id
                sess.identity_name = id_auth.identity.name
                flask_session["session"] = serialize_session(sess)
                # Attach identity to the current session
                response_object = {
                    'status': 'success',
                    'message': ""
                }
                return make_response(jsonify(response_object)), 200
            else:
                response_object = {
                    'status': 'error',
                    'message': 'Identity disabled, cannot login'
                }
                return make_response(jsonify(response_object)), 401
        else:
            response_object = {
                'status': 'error',
                'message': 'Identity not authorized'
            }
            return make_response(jsonify(response_object)), 401


# Special behavior: "authn" is a singleton, which can be None or defined with a login
# POST not implemented
view_func = AuthnAPI.as_view("authn")
bp_auth.add_url_rule(f"{app_api_base}/authn", view_func=view_func, methods=['GET', 'PUT', 'DELETE'])


@n_session()
def get_user_functions(prefix):
    """
    Elaborate a list of dictionaries, with the functions the current user can operate with from the API,
    considering ACLs for system functions
    """
    identity = g.n_session.identity
    session = g.n_session.db_session
    ahora = datetime.now()
    # Obtain system functions, and check ACLs for the identity
    functions = session.query(SystemFunction).filter(SystemFunction.name.startswith(prefix+"-")).all()
    function_uuids = {str(f.uuid): f.name for f in functions}
    function_unspecified = {f.name for f in functions}
    obj_type = session.query(ObjectType).filter(ObjectType.name == "sys-functions").first()
    # Rules stored in the database. Find the active one for each of the desired functions
    q = session.query(ACLExpression). \
        filter(and_(
                or_(ACLExpression.validity_start == None, ACLExpression.validity_start <= ahora),
                or_(ACLExpression.validity_end == None, ACLExpression.validity_end > ahora))). \
        join(ACLExpression.acl).filter(
        and_(ACL.object_type == obj_type.id, ACL.object_uuid.in_(function_uuids.keys())))
    # sql = str(q.statement.compile(dialect=postgresql.dialect()))
    rules = q.all()
    content = []
    for acl_expression in rules:
        rule = acl_expression.expression
        if rule:
            ast = string_to_ast(authr_expression, rule)
            can_execute = ast_evaluator(ast, identity)
        else:
            can_execute = False  # If no rule is found, allow execution for everybody!
        name = function_uuids[str(acl_expression.acl.object_uuid)]
        function_unspecified.discard(name)
        if can_execute:
            content.append(dict(name=name, permissions=["read"]))
        else:
            content.append(dict(name=name, permissions=[]))
    # Add the unspecified functions
    for name in function_unspecified:
        content.append(dict(name=name, permissions=["unspecified"]))

    return ResponseObject(issues=[], status=200, content=content).get_response()


bp_auth.add_url_rule(f"{app_api_base}/authn/functions/<prefix>", view_func=get_user_functions, methods=['GET'])
