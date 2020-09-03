"""
Perhaps, although functionality can be here, move authentication/authorization
to NGINX so it can filter ALL requests. An example seemingly interesting project:
https://github.com/mbreese/subauth
"""
from flask import request, abort
from functools import wraps
import firebase_admin
from firebase_admin import credentials, auth

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
