from flask import Blueprint

bp_analyses = Blueprint('analyses', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


# Analyses API

class AnalysesAPI(MethodView):
    """
    Analyses Resource
    """
    def get(self, id=None):
        print(f'GET {request.path}\nGetting analyses {id}')
        self._check_data()
        from biobarcoding.services.x import x
        response, code = x()
        return response, code


    def post(self):
        print(f'POST {request.path}\nCreating analyses')
        self._check_data()
        from biobarcoding.services.x import x
        response, code = x()
        return response, code


    def delete(self, id):
        print(f'DELETE {request.path}\nDeleting analyses {id}')
        self._check_data()
        from biobarcoding.services.x import x
        response, code = x()
        return response, code


    def _check_data(self):
        post_data = request.get_json()
        print(f'JSON data: {post_data}')


analyses = AnalysesAPI.as_view('ansis_api')
bp_analyses.add_url_rule(
    bcs_api_base + '/bos/analyses',
    view_func=analyses,
    methods=['GET','POST']
)
bp_analyses.add_url_rule(
    bcs_api_base + '/bos/analyses/<int:id>',
    view_func=analyses,
    methods=['GET','DELETE']
)
