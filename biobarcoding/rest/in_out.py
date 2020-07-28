from flask import Blueprint

bp_io = Blueprint('io', __name__)

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


file = FileAPI.as_view('file_api')

bp_io.add_url_rule(
    '/io/export',
    view_func=file,
    methods=['GET']
)
