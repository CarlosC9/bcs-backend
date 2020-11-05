"""
Perhaps, although functionality can be here, move authentication/authorization
to NGINX so it can filter ALL requests. An example seemingly interesting project:
https://github.com/mbreese/subauth
"""
import collections
import json
import sys
import traceback

import blosc
import jsonpickle
from flask import request, abort, Response, current_app, session as flask_session, g
from functools import wraps
import firebase_admin
from firebase_admin import credentials, auth
from dataclasses import dataclass, field
from typing import List, Dict, Any
import numpy as np
from multidict import MultiDict, CIMultiDict
from sqlalchemy import and_

from biobarcoding.common import generate_json
from biobarcoding.common.helpers import serialize_from_object, deserialize_to_object
from biobarcoding.db_models import DBSession
from biobarcoding.db_models.sysadmin import Authenticator, Identity, IdentityAuthenticator


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


EXEMPT_METHODS = set(['OPTIONS'])
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
                # TODO Find user with corresponding API key (pair identity-authentication method)
                pass
            elif tok_dict["auth_method"] == "local":
                # TODO Find user with USER name (pair identity-authentication method)
                pass
        elif "firebase" in tok_dict:
            # Firebase token, subauthenticator
            auth_type = session.query(Authenticator).filter(Authenticator.name == "firebase").first()
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
            # session.commit()

        iden_authenticator = session.query(IdentityAuthenticator).filter(and_(IdentityAuthenticator.identity == iden, IdentityAuthenticator.authenticator == auth_type)).first()
        if not iden_authenticator:
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

    email = request.args.get("user")
    if email:
        if tok is None:
            tok = {"email": email, "auth_method": "local"}

    if len(tok) == 0:
        abort(401, 'No authentication information provided')

    iden_authenticator = prepare_identity(tok)

    return iden_authenticator


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


class bcs_session(object):
    """
    Decorator for methods requiring: Session, and Authentication/Authorization
    @bcs_session(read_only=False)
    """
    def __init__(self, read_only=False, allowed_roles=None, disallowed_roles=[]):
        self.read_only = read_only
        self.allowed_roles = allowed_roles
        self.disallowed_roles = disallowed_roles

    def __call__(self, f):
        def wrapped_f(*args):
            sess = deserialize_session(flask_session.get("session"))
            if isinstance(sess, BCSSession):
                try:
                    g.bcs_session = sess
                    # TODO Check that identity can execute the function
                    res = f(*args)
                except:
                    res = None
                finally:
                    if not self.read_only:
                        flask_session["session"] = serialize_session(g.bcs_session)
                    # g.bcs_session = None  -- NO need for this line, "g" is reset after every request
                return res
            else:
                return sess

        return wrapped_f

