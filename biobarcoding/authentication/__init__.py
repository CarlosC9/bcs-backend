"""
Perhaps, although functionality can be here, move authentication/authorization
to NGINX so it can filter ALL requests. An example seemingly interesting project:
https://github.com/mbreese/subauth
"""
import inspect
import traceback
from dataclasses import dataclass
from datetime import datetime
from functools import wraps
from typing import Dict, Tuple, List

import blosc
import firebase_admin
from firebase_admin import credentials, auth
from flask import request, abort, Response, session as flask_session, g
from sqlalchemy import and_, or_
from sqlalchemy.orm import Session
from werkzeug.security import check_password_hash

from ..authorization import ast_evaluator, authr_expression, string_to_ast
from ..common import generate_json
from ..common.helpers import serialize_from_object, deserialize_to_object
from ..db_models import DBSession, ObjectType, DBSessionChado, DBSessionGeo
from ..db_models.sysadmin import Authenticator, Identity, IdentityAuthenticator, ACLExpression, ACL, SystemFunction, \
    Role, RoleIdentity


def initialize_firebase(app):
    cert_path = app.config['GOOGLE_APPLICATION_CREDENTIALS']
    cred = credentials.Certificate(cert_path)
    try:
        firebase_app = firebase_admin.initialize_app(cred)
    except Exception as e:
        pass


def token_required(func):
    @wraps(func)
    def decorated(*args, **kwargs):
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
        return func(*args, **kwargs)

    return decorated


"""
login -> starts session, which is serialized
logout -> ends session, 
"""


# ------------------------------


def build_json_response(obj, status=200):
    return Response(generate_json(obj),
                    mimetype="text/json",
                    status=status)


@dataclass
class AppSession:
    identity_id: int = 0
    identity_name: str = "anonymous"
    login_time: str = ""
    token_type: str = None
    token: str = None
    roles = None  # Roles requested through API Key authentication
    # Not persisted
    identity: Identity = None
    db_session: Session = None
    chado_db_session: Session = None
    postgis_db_session: Session = None


EXEMPT_METHODS = {'OPTIONS'}
NO_SESS_RESPONSE = build_json_response({"error": "No active session. Please, open one first ('PUT /api/authn')"}, 400)


def obtain_idauth_from_request() -> Tuple[IdentityAuthenticator, List[str]]:
    from ..rest import logger

    def prepare_identity(tok_dict: Dict) -> Tuple[IdentityAuthenticator, List[str]]:
        def initialize_identity_roles():
            from ..rest import tm_default_users

            # Role initialization applies to non-system users
            if str(identity.uuid) not in tm_default_users.keys():
                # if "sys-admin" role is empty (of normal identities),
                # add current identity is a normal identity, add identity to sys-admin
                sys_admin_role = session.query(Role).filter(Role.uuid == 'c79b4ff7-9576-45f6-a439-551ac23c563b').first()
                found = False
                for rol_ident in sys_admin_role.identities:
                    if str(rol_ident.identity.uuid) not in tm_default_users.keys():
                        found = True
                        break
                if not found:
                    rol_ident = RoleIdentity(identity=identity, role=sys_admin_role)
                    session.add(rol_ident)
                    # sys_admin_role.identities.append(identity)
                    # identity.roles.append(sys_admin_role)
                # if current identity does not have a role, add it to guest
                if len(identity.roles) == 0:
                    guest_role = session.query(Role).filter(Role.uuid == 'feb13f20-4223-4602-a195-a3ea14615982').first()
                    rol_ident = RoleIdentity(identity=identity, role=guest_role)
                    session.add(rol_ident)
                    # guest_role.identities.append(identity)
                    # identity.roles.append(guest_role)

        # Get Authenticator, Name and e-mail (e-mail is used as identifier in "Identity")
        identity = None
        authenticator = None
        ident_auth = None
        name = None
        email = None
        roles = None  # Not specified (not empty)
        session = DBSession()
        if "auth_method" in tok_dict:
            if tok_dict["auth_method"] == "local-api-key":
                authenticator = session.query(Authenticator).filter(Authenticator.name == "local-api-key").first()
                identity = session.query(Identity).filter(Identity.name == tok_dict["user"]).first()
                if identity:
                    ident_auth = session.query(IdentityAuthenticator).\
                        filter(and_(IdentityAuthenticator.identity_id == identity.id,
                                    IdentityAuthenticator.authenticator_id == authenticator.id)).first()
                    if ident_auth:
                        found = False
                        for d in ident_auth.authenticator_info:
                            if check_password_hash(d["hash"], tok_dict["api_key"]):
                                found = True
                                # Check if dates are valid
                                if "valid_from" in d:
                                    if d["valid_from"] is not None:
                                        d["valid_from"] = datetime.fromisoformat(d["valid_from"])
                                else:
                                    d["valid_from"] = None
                                if "valid_until" in d:
                                    if d["valid_until"] is not None:
                                        d["valid_until"] = datetime.fromisoformat(d["valid_until"])
                                else:
                                    d["valid_until"] = None
                                if (d["valid_from"] is None or (d["valid_from"] is not None and d["valid_from"] <= datetime.now())) and (d["valid_until"] is None or (d["valid_until"] is not None and d["valid_until"] >= datetime.now())):
                                    roles = d["roles"]
                                    # TODO - Intersect with current roles of Identity
                                    email = identity.email
                                    name = identity.name
                                else:
                                    logger.warning("Specified API key is not in the range of valid dates")
                                break
                    if roles is None:
                        logger.warning("Specified API key is not valid")
                    else:
                        identity = None

            elif tok_dict["auth_method"] == "local":
                authenticator = session.query(Authenticator).filter(Authenticator.name == "local").first()
                name = tok_dict["user"]
                # TODO Check the network the request is coming from.
                #  It must be the same network of the backend
                request_ip = request.remote_addr
        elif "firebase" in tok_dict:
            # Firebase token, subauthenticator
            authenticator = session.query(Authenticator).filter(Authenticator.name == "firebase").first()
            if tok["firebase"]["sign_in_provider"] == "anonymous":
                name = "_anonymous"
            elif tok["firebase"]["sign_in_provider"] == "password":
                email = tok_dict["email"]
            else:
                name = tok_dict["name"]
                email = tok_dict["email"]

        # Find or create Identity
        if not identity:
            if name:
                identity = session.query(Identity).filter(Identity.name == name).first()
        if not identity:
            identity = session.query(Identity).filter(Identity.email == email).first()
        if not identity:
            identity = Identity()
            identity.can_login = True
            identity.email = email
            identity.name = name
            session.add(identity)

        # Find or create IdentityAuthenticator
        if not ident_auth:
            ident_auth = session.query(IdentityAuthenticator).filter(
                and_(IdentityAuthenticator.identity == identity, IdentityAuthenticator.authenticator == authenticator)).first()
        if not ident_auth and authenticator.name == "firebase":  # Create (others must exist previously)
            ident_auth = IdentityAuthenticator()
            ident_auth.email = email
            ident_auth.name = name
            ident_auth.authenticator_info = tok_dict
            ident_auth.identity = identity
            ident_auth.authenticator = authenticator
            session.add(ident_auth)

        if ident_auth:
            ident_auth.last_login_time = datetime.utcnow()

        # Roles initialization
        initialize_identity_roles()

        session.commit()

        return ident_auth, roles

    # ------------------------------------------------------------------------------------------------------------------
    # Get the token from the request
    tok = {}
    if 'Authorization' in request.headers:  # Firebase authentication
        auth_token = request.headers['Authorization'].split(" ")[1]
        try:
            tok = auth.verify_id_token(auth_token)
        except Exception as e:
            logger.debug(e)
            traceback.print_exc()
            abort(401, 'The session token is not valid or has expired')
    elif 'X-API-Key' in request.headers:  # "API-KEY" authentication
        user = request.args.get("user")
        tok = dict(api_key=request.headers["X-API-Key"], auth_method="local-api-key", user=user)
    elif "user" in request.args:  # "local" authentication
        user = request.args.get("user")
        tok = dict(user=user, auth_method="local")

    if len(tok) == 0:
        abort(401, 'No authentication information provided')

    return prepare_identity(tok)


def serialize_session(state: AppSession):
    tmp = serialize_from_object(state)
    tmp = blosc.compress(bytearray(tmp, "utf-8"), cname="zlib", typesize=8)
    return tmp


def deserialize_session(s, return_error_response_if_none=True) -> AppSession:
    def deserialize_session_sub(st: str) -> AppSession:
        if isinstance(st, bytes):
            st = blosc.decompress(st).decode("utf-8")
        if isinstance(st, str):
            state = deserialize_to_object(st)
        else:
            raise Exception(f"Serialized session must be a string: currently is of type {type(st)}")

        return state

    if s:
        try:
            sess = deserialize_session_sub(s)
        except Exception as e:
            traceback.print_exc()
            sess = None
    else:
        sess = None

    if not sess and return_error_response_if_none:
        return NO_SESS_RESPONSE
    else:
        return sess


class n_session(object):
    """
    Decorator for RESTful methods requiring: Session, and Authentication/Authorization
    @n_session(read_only=False)
    """

    def __init__(self, read_only=False, authr=None):
        """
        The authentication "can execute" rule can be:
          - None: assume the function name, using "reflection"
          - Simple string: a literal with a function name regarding permissions
          - String following an "authr_expression" Authorization Rule syntax
        :param read_only: True if it is a read only function, regarding database
        :param authr: Authentication rule
        """
        self.read_only = read_only
        self.authr = authr

    @staticmethod
    def is_function_name(s: str):
        t = s.split(" ")
        if len(t) < 4:
            valid = True
            for p in t:
                if not p.isidentifier():
                    valid = False
        else:
            valid = False

        return valid

    def __call__(self, f):
        def wrapped_f(*args, **kwargs):

            def identity_can_impersonate(iden: Identity):
                return iden.name in ["celery_user", "test_user"]

            sess = deserialize_session(flask_session.get("session"))
            if isinstance(sess, AppSession):
                try:
                    g.n_session = sess
                    db_session = DBSession()
                    g.n_session.db_session = db_session
                    ident = db_session.query(Identity).get(g.n_session.identity_id)
                    # Impersonate if:
                    # - Impersonate header is set
                    # - Current user can impersonate (e.g. celery_user)
                    # - Operation is read-only
                    # Example:
                    #   curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=celery_user"
                    #   curl --cookie app-cookies.txt "$API_BASE_URL/geo/layers/" -H "Impersonated-id: 30"
                    if 'Impersonated-id' in request.headers and identity_can_impersonate(ident) and self.read_only:
                        g.n_session.original_identity_id = g.n_session.identity_id
                        g.n_session.identity_id = request.headers['Impersonated-id']
                        ident = db_session.query(Identity).get(g.n_session.identity_id)
                    g.n_session.identity = ident
                    chado_db_session = DBSessionChado()  # Chado test
                    g.n_session.chado_db_session = chado_db_session  # Chado test
                    postgis_db_session = DBSessionGeo()
                    g.n_session.postgis_db_session = postgis_db_session
                    try:  # Protect "db_session"
                        # Get execution rule
                        if self.authr is None or n_session.is_function_name(self.authr):
                            if self.authr is None:
                                # (be careful of refactorizations)
                                f_code_name = f"{inspect.getmodule(f).__name__}.{f.__name__}"
                            else:
                                f_code_name = self.authr
                            # Search authorization database "f_code_name"
                            ahora = datetime.now()
                            function = db_session.query(SystemFunction).filter(
                                SystemFunction.name == f_code_name).first()
                            if function:
                                obj_type = db_session.query(ObjectType).filter(
                                    ObjectType.name == "sys-functions").first()
                                # Rule stored in the database. Find the active one for the desired function
                                # TODO if there is no ACLExpression, search ACL elements
                                #  (or maybe generate an expression from the ACL elements)
                                rule = db_session.query(ACLExpression.expression). \
                                    filter(and_(
                                    or_(ACLExpression.validity_start == None, ACLExpression.validity_start <= ahora),
                                    or_(ACLExpression.validity_end == None, ACLExpression.validity_end > ahora))). \
                                    join(ACLExpression.acl).filter(
                                    and_(ACL.object_type == obj_type.id, ACL.object_uuid == function.uuid)).first()
                            else:
                                rule = None  # If the function name cannot be found, rule = None -> so "can execute"
                        else:
                            rule = self.authr  # Rule specified literally, needs to be parsed
                        # Check execution permission
                        if rule:
                            ast = string_to_ast(authr_expression, rule)
                            can_execute = ast_evaluator(ast, ident)
                        else:
                            can_execute = True  # If no rule is found, allow execution for everybody!
                        # Execute
                        if can_execute:
                            g.commit_after = True
                            res = f(*args, **kwargs)
                            if g.commit_after:
                                db_session.commit()
                                chado_db_session.commit()  # Chado test
                                postgis_db_session.commit()
                        else:
                            raise Exception(
                                f"Current user ({sess.identity_name}) cannot execute function {f_code_name}")
                    except:
                        traceback.print_exc()
                        db_session.rollback()
                        chado_db_session.rollback()  # Chado test
                        postgis_db_session.rollback()
                        raise
                    finally:
                        DBSessionGeo.remove()
                        DBSession.remove()
                        DBSessionChado.remove()  # Chado test
                except:
                    res = None
                finally:
                    # Restore impersonating identity
                    if hasattr(g.n_session, "original_identity_id") and g.n_session.original_identity_id:
                        g.n_session.identity_id = g.n_session.original_identity_id
                        g.n_session.original_identity_id = None

                    if not self.read_only:
                        g.n_session.identity = None
                        g.n_session.db_session = None
                        g.n_session.chado_db_session = None  # Chado test
                        g.n_session.postgis_db_session = None
                        flask_session["session"] = serialize_session(g.n_session)
                    # g.n_session = None  -- NO need for this line, "g" is reset after every request
                return res
            else:
                return sess

        return wrapped_f
