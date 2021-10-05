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
from typing import Dict

import blosc
import firebase_admin
from firebase_admin import credentials, auth
from flask import request, abort, Response, session as flask_session, g
from sqlalchemy import and_, or_
from sqlalchemy.orm import Session

from ..authorization import ast_evaluator, authr_expression, string_to_ast
from ..common import generate_json
from ..common.helpers import serialize_from_object, deserialize_to_object
from ..db_models import DBSession, ObjectType, DBSessionChado, DBSessionGeo
from ..db_models.sysadmin import Authenticator, Identity, IdentityAuthenticator, ACLExpression, ACL, \
    SystemFunction


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
class BCSSession:
    identity_id: int = 0
    identity_name: str = "anonymous"
    login_time: str = ""
    token_type: str = None
    token: str = None
    # Not persisted
    identity: Identity = None
    db_session: Session = None
    chado_db_session: Session = None
    postgis_db_session: Session = None


EXEMPT_METHODS = {'OPTIONS'}
NO_SESS_RESPONSE = build_json_response({"error": "No active session. Please, open one first ('PUT /api/authn')"}, 400)


# def login_required(func):
#     @wraps(func)
#     def decorated_view(*args, **kwargs):
#         if request.method in EXEMPT_METHODS:
#             return func(*args, **kwargs)
#         elif current_app.config.get('LOGIN_DISABLED'):
#             return func(*args, **kwargs)
#         elif current_session is None:
#             return NO_SESS_RESPONSE
#         return func(*args, **kwargs)
#     return decorated_view


def obtain_idauth_from_request() -> str:
    def prepare_identity(tok_dict: Dict) -> IdentityAuthenticator:
        name = None
        email = None
        auth_type = None
        session = DBSession()
        if "auth_method" in tok_dict:
            if tok_dict["auth_method"] == "local-api-key":
                auth_type = session.query(Authenticator).filter(Authenticator.name == "local-api-key").first()
            elif tok_dict["auth_method"] == "local":
                auth_type = session.query(Authenticator).filter(Authenticator.name == "local").first()
                name = tok_dict["user"]
        elif "firebase" in tok_dict:
            # Firebase token, subauthenticator
            auth_type = session.query(Authenticator).filter(Authenticator.name == "firebase").first()
            if tok["firebase"]["sign_in_provider"] == "anonymous":
                name = "_anonymous"
            else:
                name = tok_dict["name"]
                email = tok_dict["email"]
        iden = session.query(Identity).filter(Identity.name == name).first()
        if not iden:
            iden = session.query(Identity).filter(Identity.email == email).first()
        if not iden:
            iden = Identity()
            iden.email = email
            iden.name = name
            session.add(iden)

        iden_authenticator = session.query(IdentityAuthenticator).filter(
            and_(IdentityAuthenticator.identity == iden, IdentityAuthenticator.authenticator == auth_type)).first()
        if not iden_authenticator and auth_type.name == "firebase":  # Create (others must exist previously)
            iden_authenticator = IdentityAuthenticator()
            iden_authenticator.email = email
            iden_authenticator.name = name
            iden_authenticator.authenticator_info = tok_dict
            iden_authenticator.identity = iden
            iden_authenticator.authenticator = auth_type
            session.add(iden_authenticator)
        session.commit()
        return iden_authenticator

    tok = None
    if 'Authorization' in request.headers:
        auth_token = request.headers['Authorization'].split(" ")[1]
        try:
            tok = auth.verify_id_token(auth_token)
            print(tok)
        except Exception as e:
            print(e)
            abort(401, 'The session token is not valid or has expired')
    elif 'X-API-Key' in request.headers:
        # api_key
        tok = dict(api_key=request.headers["X-API-Key"], auth_method="local-api-key")

    user = request.args.get("user")
    if user:
        if tok is None:
            tok = {"user": user, "auth_method": "local"}

    if len(tok) == 0:
        abort(401, 'No authentication information provided')

    return prepare_identity(tok)


def serialize_session(state: BCSSession):
    tmp = serialize_from_object(state)
    tmp = blosc.compress(bytearray(tmp, "utf-8"), cname="zlib", typesize=8)
    return tmp


def deserialize_session(s, return_error_response_if_none=True) -> BCSSession:
    def deserialize_session_sub(st: str) -> BCSSession:
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
    @bcs_session(read_only=False)
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
            sess = deserialize_session(flask_session.get("session"))
            if isinstance(sess, BCSSession):
                try:
                    g.n_session = sess
                    db_session = DBSession()
                    g.n_session.db_session = db_session
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
                                # TODO if there is no ACLExpression, search ACL elements (or maybe generate an expression from the ACL elements)
                                rule = db_session.query(ACLExpression.expression). \
                                    filter(and_(
                                    or_(ACLExpression.validity_start is None, ACLExpression.validity_start <= ahora),
                                    or_(ACLExpression.validity_end is None, ACLExpression.validity_end > ahora))). \
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
                            res = f(*args, **kwargs)
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
                    if not self.read_only:
                        g.n_session.identity = None
                        g.n_session.db_session = None
                        g.n_session.chado_db_session = None  # Chado test
                        g.n_session.postgis_db_session = None
                        flask_session["session"] = serialize_session(g.n_session)
                    # g.bcs_session = None  -- NO need for this line, "g" is reset after every request
                return res
            else:
                return sess

        return wrapped_f
