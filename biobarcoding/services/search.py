from flask import request, make_response, jsonify
from flask.views import MethodView

class SearchAPI(MethodView):
    """
    Search Resource
    """
    def get(self, type=None):
        msg = f'GET {request.path}\nGetting search {type}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def _check_data(self):
        post_data = request.get_json()
        print(f'JSON data: {post_data}')
