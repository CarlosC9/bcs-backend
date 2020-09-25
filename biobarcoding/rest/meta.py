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
        print(f'GET {request.path}\nGetting taxonomy {id}')
        self._check_data()
        from biobarcoding.services.taxonomies import read_taxonomy
        resp = read_taxonomy()
        responseObject = {
            'status': 'success',
            'message': resp
        }
        return make_response(jsonify(responseObject)), 200


    @token_required
    def post(self):
        print(f'POST {request.path}\nCreating taxonomy')
        self._check_data()
        if not self.name:
            return make_response('Parameters missing: "name" is required.', 400)
        from biobarcoding.services.taxonomies import create_taxonomy
        resp = create_taxonomy(self.name)
        responseObject = {
            'status': 'success',
            'message': resp
        }
        return make_response(jsonify(responseObject)), 200


    @token_required
    def delete(self, id):
        print(f'DELETE {request.path}\nDeleting taxonomy {id}')
        self._check_data()
        from biobarcoding.services.taxonomies import delete_taxonomy
        resp = delete_taxonomy(id)
        responseObject = {
        'status': 'success',
        'message': resp
        }
        return make_response(jsonify(responseObject)), 200


    def _check_data(self):
        post_data = request.get_json()
        self.name = None
        if post_data and 'name' in post_data:
            self.name = post_data['name']
        print(f'JSON data: {post_data}')


taxon = TaxonomyAPI.as_view('taxon_api')
bp_meta.add_url_rule(
    bcs_api_base + '/bo/taxonomy',
    view_func=taxon,
    methods=['GET','POST']
)
bp_meta.add_url_rule(
    bcs_api_base + '/bo/taxonomy/<int:id>',
    view_func=taxon,
    methods=['GET','DELETE']
)


# Organism API

class OrganismAPI(MethodView):
    """
    Organism Resource
    """
    @token_required
    def get(self, id=None):
        print(f'GET {request.path}\nGetting organism {id}')
        self._check_data()
        from biobarcoding.services.organisms import read_organism
        resp = read_organism()
        responseObject = {
            'status': 'success',
            'message': resp
        }
        return make_response(jsonify(responseObject)), 200


    @token_required
    def post(self):
        print(f'POST {request.path}\nCreating organism')
        self._check_data()
        if not ('genus' in request.json and 'species' in request.json):
            responseObject = {
                'status': 'error',
                'message': 'Parameters missing: genus and species are required.'
            }
            return make_response(jsonify(responseObject)), 400
        from biobarcoding.services.organisms import create_organism
        resp = create_organism(request.json['genus'],request.json['species'])
        responseObject = {
            'status': 'success',
            'message': resp
        }
        return make_response(jsonify(responseObject)), 200


    @token_required
    def delete(self, id):
        print(f'DELETE {request.path}\nDeleting organism {id}')
        self._check_data()
        from biobarcoding.services.organisms import delete_organism
        resp = delete_organism(id)
        responseObject = {
        'status': 'success',
        'message': resp
        }
        return make_response(jsonify(responseObject)), 200


    def _check_data(self):
        post_data = request.get_json()
        print(f'JSON data: {post_data}')


organism = OrganismAPI.as_view('organism_api')
bp_meta.add_url_rule(
    bcs_api_base + '/bo/organism',
    view_func=organism,
    methods=['GET','POST']
)
bp_meta.add_url_rule(
    bcs_api_base + '/bo/organism/<int:id>',
    view_func=organism,
    methods=['GET','DELETE']
)


# Analysis API

class AnalysisAPI(MethodView):
    """
    Analysis Resource
    """
    def get(self, id=None):
        print(f'GET {request.path}\nGetting analysis {id}')
        self._check_data()
        from biobarcoding.services.analyses import read_analysis
        resp = read_analysis()
        responseObject = {
            'status': 'success',
            'message': resp
        }
        return make_response(jsonify(responseObject)), 200


    def post(self):
        print(f'POST {request.path}\nCreating analysis')
        self._check_data()
        if not ('program' in request.json and 'programversion' in request.json):
            responseObject = {
                'status': 'error',
                'message': 'Parameters missing: program and programversion are required.'
            }
            return make_response(jsonify(responseObject)), 400
        from biobarcoding.services.analyses import create_analysis
        resp = create_analysis(request.json['program'],request.json['programversion'])
        responseObject = {
            'status': 'success',
            'message': resp
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        print(f'DELETE {request.path}\nDeleting analysis {id}')
        self._check_data()
        from biobarcoding.services.analyses import delete_analysis
        resp = delete_analysis(id)
        responseObject = {
        'status': 'success',
        'message': resp
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
        print(f'Getting ontology {id}')
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': resp
        }
        return make_response(jsonify(responseObject)), 200


    def post(self):
        print(f'Creating ontology')
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': resp
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        print(f'Deleting ontology {id}')
        self._check_data()

        responseObject = {
        'status': 'success',
        'message': resp
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
