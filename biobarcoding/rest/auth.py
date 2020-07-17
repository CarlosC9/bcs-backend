from flask import Blueprint
from biobarcoding.rest import bcs_api_base

bp_auth = Blueprint('auth', __name__)

from biobarcoding.services import AuthAPI

auth = AuthAPI.as_view('auth_api')
bp_auth.add_url_rule(
    bcs_api_base + '/auth/login',
    view_func=auth,
    methods=['POST']
)
bp_auth.add_url_rule(
    bcs_api_base + '/auth/logout',
    view_func=auth,
    methods=['POST']
)
