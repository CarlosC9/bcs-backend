from flask import Blueprint

bp_organisms = Blueprint('organisms', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


# Organisms API

class OrganismsAPI(MethodView):
    """
    Organisms Resource
    """
    def get(self, id=None):
        print(f'GET {request.path}\nGetting organisms {id}')
        self._check_data()
        from biobarcoding.services.x import x
        response, code = x()
        return response, code


    def post(self):
        print(f'POST {request.path}\nCreating organisms')
        self._check_data()
        from biobarcoding.services.x import x
        response, code = x()
        return response, code


    def delete(self, id):
        print(f'DELETE {request.path}\nDeleting organisms {id}')
        self._check_data()
        from biobarcoding.services.x import x
        response, code = x()
        return response, code


    def _check_data(self):
        post_data = request.get_json()
        print(f'JSON data: {post_data}')


organisms = OrganismsAPI.as_view('taxon_api')
bp_organisms.add_url_rule(
    bcs_api_base + '/bos/organisms',
    view_func=organisms,
    methods=['GET','POST']
)
bp_organisms.add_url_rule(
    bcs_api_base + '/bos/organisms/<int:id>',
    view_func=organisms,
    methods=['GET','DELETE']
)
