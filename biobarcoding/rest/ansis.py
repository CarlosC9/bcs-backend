from flask import Blueprint

bp_ansis = Blueprint('analysis', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

class AnalysisAPI(MethodView):
    """
    Analysis Resource
    """
    def get(self, id=None):
        msg = f'GET {request.path}\nGetting analysis {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def post(self):
        msg = f'POST {request.path}\nCreating analysis'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting analysis {id}'
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


ansis = AnalysisAPI.as_view('ansis_api')
bp_ansis.add_url_rule(
    '/bo/analysis',
    view_func=ansis,
    methods=['GET','POST']
)
bp_ansis.add_url_rule(
    '/bo/analysis/<int:id>',
    view_func=ansis,
    methods=['GET','DELETE']
)
