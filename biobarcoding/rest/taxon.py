from flask import Blueprint

bp_taxon = Blueprint('taxonomy', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

class TaxonomyAPI(MethodView):
    """
    Taxonomy Resource
    """
    @token_required
    def get(self, id=None):
        msg = f'GET {request.path}\nGetting taxonomy {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200

    @token_required
    def post(self):
        msg = f'POST {request.path}\nCreating taxonomy'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200

    @token_required
    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting taxonomy {id}'
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


taxon = TaxonomyAPI.as_view('taxon_api')
bp_taxon.add_url_rule(
    '/bo/organism',
    view_func=taxon,
    methods=['GET','POST']
)
bp_taxon.add_url_rule(
    '/bo/organism/<int:id>',
    view_func=taxon,
    methods=['GET','DELETE']
)
