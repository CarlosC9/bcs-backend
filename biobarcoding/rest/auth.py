from flask import Blueprint
from biobarcoding.rest import bcs_api_base

bp_auth = Blueprint('auth', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

class AuthAPI(MethodView):
    """
    Authentication Resource
    """
    def post(self):
        msg = f'GET {request.path}\nGetting analysis {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def _check_data(self):
        if request.path == '/auth/login':
            msg = 'Successfully logged in.'
        elif request.path == '/auth/logout':
            msg = 'Successfully logged in.'
        else:
            msg = 'Warning: Route may be wrong.'

        post_data = request.get_json()
        print(f'{msg}\nJSON data: {post_data}')


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
