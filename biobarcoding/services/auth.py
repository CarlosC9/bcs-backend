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
