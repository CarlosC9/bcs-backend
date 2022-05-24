import secrets
import uuid
from datetime import datetime

from firebase_admin import auth
from flask import Blueprint, abort, request, make_response, jsonify, session as flask_session, g
from flask.views import MethodView
from sqlalchemy import and_, or_
from werkzeug.security import generate_password_hash

from ..authentication import serialize_session, AppSession, deserialize_session, obtain_idauth_from_request, n_session
from . import app_api_base, ResponseObject, logger, Issue, IType, register_api
from ..authorization import string_to_ast, authr_expression, ast_evaluator
from ..common.helpers import obj_to_json
from ..db_models import ObjectType
from ..db_models.sysadmin import SystemFunction, ACLExpression, ACL, Identity, Authenticator, IdentityAuthenticator, \
    Role

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


bp_api_key = Blueprint('api_key', __name__)


class ApiKeyAPI(MethodView):
    """
    To manage API keys of users

    LOGIN (entries assume the current user):
    export API_BASE_URL=http://localhost:5000/api
    curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"

    LIST USERS WITH API KEYS:
    curl --cookie app-cookies.txt "$API_BASE_URL/api_keys/"
    GET API KEYS FOR A USER:
    curl --cookie app-cookies.txt "$API_BASE_URL/api_keys/<identity_id>"
    GET API KEYS FOR CURRENT USER:
    curl --cookie app-cookies.txt "$API_BASE_URL/api_keys/0"

    CREATE API KEY for CURRENT USER:
    curl --cookie app-cookies.txt -X POST "$API_BASE_URL/api_keys/" -H "Content-Type: application/json" -d '{"roles": ["read-molecular-api"]}'
    TODO - CREATE API Key with a period of validity
    DELETE API KEY for CURRENT USER:
    * First GET API KEYS for CURRENT USER, take note of key_idx field of the Key to delete
    curl --cookie app-cookies.txt -X DELETE "$API_BASE_URL/api_keys/<key_idx>"

    EXAMPLE
    LOGIN with API KEY:
    curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user" -H "X-API-Key: 8217b03ac2f34cfabd1388d28d420387"

    """
    @n_session()
    def get(self, identity_id=None):
        """
        Obtain the list of API keys of a user, or the users having API key authentication enabled

        :param identity_id:
        :return:
        """
        session = g.n_session.db_session
        authenticator = session.query(Authenticator).filter(Authenticator.name == "local-api-key").first()
        r = ResponseObject()
        if identity_id is not None:
            if identity_id == 0:  # Current user
                identity_id = g.n_session.identity.id
            # Obtain the list of API keys of a user
            identity = session.query(Identity).filter(Identity.id == identity_id).first()
            if not identity:
                r.issues.append(Issue(IType.ERROR, f"Identity {identity_id} not found"))
            else:
                # Read list of API Keys from the relation with "API Key Authenticator", do not show the hash (does not make sense)
                identity_authenticator = session.query(IdentityAuthenticator).\
                    filter(and_(IdentityAuthenticator.identity_id == identity_id,
                                IdentityAuthenticator.authenticator_id == authenticator.id)).first()
                if not identity_authenticator:
                    r.issues.append(Issue(IType.WARNING, f"Identity {identity_id} does not have API Key authentication enabled"))
                else:
                    r.content = [dict(roles=d['roles'], key_idx=d["key_idx"],
                                      valid_from=d["valid_from"], valid_until=d['valid_until'])
                                 for d in identity_authenticator.authenticator_info]
        else:
            # List all identities with API keys
            identities = session.query(IdentityAuthenticator).\
                filter(IdentityAuthenticator.authenticator_id == authenticator.id).all()
            r.content = identities
        return r.get_response()

    @n_session()
    def post(self):
        """
        Add an API key for the current user

        :return:
        """
        session = g.n_session.db_session
        authenticator = session.query(Authenticator).filter(Authenticator.name == "local-api-key").first()
        r = ResponseObject()
        t = request.json  # roles, valid_from, valid_until
        # Check that the specified roles exist at the moment of creation
        if 'roles' not in t:
            r.issues.append(Issue(IType.ERROR, "'roles' are not specified"))
            r.status = 400
        else:
            roles = session.query(Role).filter(Role.name.in_(t['roles'])).all()
            if len(roles) != len(t['roles']):
                # Find which roles do not exist and show them
                lst = []
                _ = [role.name for role in roles]
                for role in t["roles"]:
                    if role not in _:
                        lst.append(role)
                r.issues.append(Issue(IType.ERROR, f"Roles: {', '.join(lst)} do not exist"))
                r.status = 400
            else:
                if 'valid_from' not in t:
                    t['valid_from'] = datetime.now().isoformat()
                if 'valid_until' not in t:
                    t['valid_until'] = None

                identity_id = g.n_session.identity.id
                identity_authenticator = session.query(IdentityAuthenticator).\
                    filter(and_(IdentityAuthenticator.identity_id == identity_id,
                                IdentityAuthenticator.authenticator_id == authenticator.id)).first()
                if not identity_authenticator:
                    identity_authenticator = IdentityAuthenticator()
                    identity_authenticator.identity_id = identity_id
                    identity_authenticator.authenticator_id = authenticator.id
                    identity_authenticator.authenticator_info = []
                    session.add(identity_authenticator)
                # Generate key_idx (to address keys pertaining to an identity)
                used_idxs = set([d['key_idx'] for d in identity_authenticator.authenticator_info])
                key_idx = 1
                while key_idx in used_idxs:
                    key_idx += 1
                # Generate API Key
                k = uuid.uuid4().hex
                # Fields needed to login (using "i-Bond" client)
                r.content = dict(api_key=k, name=g.n_session.identity.name)
                # Store hash
                h = generate_password_hash(k)
                d = dict(key_idx=key_idx, hash=h, roles=t['roles'],
                         valid_from=t['valid_from'], valid_until=t['valid_until'])
                lst = identity_authenticator.authenticator_info.copy()
                lst.append(d)
                identity_authenticator.authenticator_info = lst
                d2 = d.copy()
                del d2['hash']
                del d2['key_idx']
                r.content.update(d2)
        return r.get_response()

    @n_session()
    def delete(self, key_idx):
        """
        Delete one of the API key's of the current user

        :param key_idx:
        :return:
        """
        session = g.n_session.db_session
        authenticator = session.query(Authenticator).filter(Authenticator.name == "local-api-key").first()
        r = ResponseObject()
        identity_id = g.n_session.identity.id
        identity_authenticator = session.query(IdentityAuthenticator).\
            filter(and_(IdentityAuthenticator.identity_id == identity_id,
                        IdentityAuthenticator.authenticator_id == authenticator.id)).first()
        if not identity_authenticator:
            r.issues.append(Issue(IType.ERROR, "No API keys found"))
            r.status = 404
        else:
            tmp = len(identity_authenticator.authenticator_info)
            identity_authenticator.authenticator_info = \
                [d for d in identity_authenticator.authenticator_info if d['key_idx'] != key_idx]
            if len(identity_authenticator.authenticator_info) == tmp:
                r.issues.append(Issue(IType.ERROR, f"API key with index {key_idx} not found"))
                r.status = 404
            else:
                r.status = 204
                if len(identity_authenticator.authenticator_info) == 0:
                    session.delete(identity_authenticator)

        return r.get_response()


url = f"{app_api_base}/api_keys/"
view_func = ApiKeyAPI.as_view("api_key")
bp_api_key.add_url_rule(url, defaults=dict(identity_id=None), view_func=view_func, methods=['GET'])
bp_api_key.add_url_rule(url, view_func=view_func, methods=['POST'])
bp_api_key.add_url_rule(f'{url}<int:identity_id>', view_func=view_func, methods=['GET'])
bp_api_key.add_url_rule(f'{url}<int:key_idx>', view_func=view_func, methods=['DELETE'])


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
        id_auth, roles = obtain_idauth_from_request()  # MAIN !!
        # If the identity is Active, continue;
        # If not, return an error
        if id_auth:
            # logger.debug(f"AUTH: {obj_to_json(id_auth)}")
            # logger.debug(f"IDENTITY: {obj_to_json(id_auth.identity)}")
            if not id_auth.identity.deactivation_time:  # and id_auth.identity.can_login
                sess = AppSession()
                sess.identity_id = id_auth.identity.id
                sess.identity_name = id_auth.identity.name
                sess.roles = roles
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
