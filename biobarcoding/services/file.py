from flask import request, make_response, jsonify
from flask.views import MethodView

class FileAPI(MethodView):
    """
    File Resource
    """
    def get(self):
        msg = f'GET {request.path}\nDownloading file.'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200
