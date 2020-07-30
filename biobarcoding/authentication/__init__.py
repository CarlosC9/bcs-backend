"""
Perhaps, although functionality can be here, move authentication/authorization
to NGINX so it can filter ALL requests. An example seemingly interesting project:
https://github.com/mbreese/subauth
"""
from flask import request, abort
from functools import wraps
import os
import firebase_admin
from firebase_admin import credentials, auth

from biobarcoding import flask_app
# cert_path = flask_app.config['GOOGLE_APPLICATION_CREDENTIALS']
# cert_path = '/home/acurbelo/.local/share/bcs-backend/firebase-key.json'
cert_path = '../firebase-key.json'
cred = credentials.Certificate(cert_path)
firebase_app = firebase_admin.initialize_app(cred)


def token_required(func):
    @wraps(func)
    def decorated(*args, **kwargs):
        auth_token = None
        if 'Authorization' in request.headers:
            auth_token = request.headers['Authorization'].split(" ")[1]
        if not auth_token:
            abort(401, 'The session token is missing')
        try:
            print(auth.verify_id_token(auth_token))
            # current_user = User.query.filter(User.id == user_id).first()
        except Exception as e:
            print(e)
            abort(401, 'The session token is not valid or has expired')
        return func(*args, **kwargs)
    return decorated

def admin_token_required(func):
    @wraps(func)
    def decorated(*args, **kwargs):
        auth_token = None
        if 'Authorization' in request.headers:
            auth_token = request.headers['Authorization'].split(" ")[1]
        if not auth_token:
            abort(401, 'The session token is missing')
        try:
            if not WhitelistToken.check_whitelist(auth_token):
                abort(401, 'The session token is not valid')
            user_id = User.decode_auth_token(auth_token)
            current_user = User.query.filter(User.id == user_id).first()
            admin = current_user.admin
            if not admin:
                abort(401, 'The session token does not have administrator permission')
        except Exception as e:
            abort(401, 'The session token is not valid or has expired')
        return func(*args, **kwargs)
    return decorated
