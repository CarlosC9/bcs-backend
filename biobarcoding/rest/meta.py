from flask import Blueprint

bp_meta = Blueprint('metadata', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base
from biobarcoding.authentication import *


# Taxonomy API

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
bp_meta.add_url_rule(
    bcs_api_base + '/bo/organism',
    view_func=taxon,
    methods=['GET','POST']
)
bp_meta.add_url_rule(
    bcs_api_base + '/bo/organism/<int:id>',
    view_func=taxon,
    methods=['GET','DELETE']
)


# Analysis API

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
bp_meta.add_url_rule(
    bcs_api_base + '/bo/analysis',
    view_func=ansis,
    methods=['GET','POST']
)
bp_meta.add_url_rule(
    bcs_api_base + '/bo/analysis/<int:id>',
    view_func=ansis,
    methods=['GET','DELETE']
)


# Ontology API

class OntologyAPI(MethodView):
    """
    Ontology Resource
    """
    def get(self, id=None):
        msg = f'Getting ontology {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def post(self):
        msg = f'Creating ontology'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        msg = f'Deleting ontology {id}'
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


onto = OntologyAPI.as_view('onto_api')
bp_meta.add_url_rule(
    bcs_api_base + '/bo/ontology',
    view_func=onto,
    methods=['GET','POST']
)
bp_meta.add_url_rule(
    bcs_api_base + '/bo/ontology/<int:id>',
    view_func=onto,
    methods=['GET','DELETE']
)
