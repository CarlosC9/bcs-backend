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

from biobarcoding.common import generate_json
from biobarcoding.common.helpers import serialize_from_object, deserialize_to_object


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


def obtain_identity_from_request() -> str:
    auth_token = None
    if 'Authorization' in request.headers:
        auth_token = request.headers['Authorization'].split(" ")[1]
    elif 'X-API-Key' in request.headers:
        # api_key
        abort(405, 'The X-API-Key authentication is not available.')
    user = request.args.get("user")
    if user:
        return user
    if not auth_token:
        abort(401, 'The session token is missing')
    try:
        print(auth.verify_id_token(auth_token))
        return auth_token
    except Exception as e:
        print(e)
        abort(401, 'The session token is not valid or has expired')


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
                    g.current_session = sess
                    # TODO Check that identity can execute the function
                    res = f(*args)
                except:
                    res = None
                finally:
                    if not self.read_only:
                        flask_session["session"] = serialize_session(g.current_session)
                    # g.current_session = None  -- NO need for this line, "g" is reset after every request
                return res
            else:
                return sess

        return wrapped_f

