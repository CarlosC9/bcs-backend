from flask import Blueprint

bp_phylotrees = Blueprint('phylo', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


class PhyloAPI(MethodView):
    """
    Phylo Resource
    """
    def get(self, id=None):
        msg = f'GET {request.path}\nGetting tree {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def post(self):
        msg = f'POST {request.path}\nCreating tree'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def put(self, id):
        msg = f'PUT {request.path}\nCreating tree {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting tree {id}'
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


class PhyloFeatAPI(MethodView):
    """
    Phylogenetic Tree Feature Resource
    """
    def get(self, id=None):
        msg = f'GET {request.path}\nGetting comment {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def post(self):
        msg = f'POST {request.path}\nCreating comment'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def put(self, id):
        msg = f'PUT {request.path}\nCreating comment {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting comment {id}'
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


phylo = PhyloAPI.as_view('phylo_api')
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bo/phylo/',
    view_func=phylo,
    methods=['GET','POST']
)
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bo/phylo/<int:phylo_id>',
    view_func=phylo,
    methods=['GET','PUT','DELETE']
)

phylo_feat = PhyloFeatAPI.as_view('phylo_feat_api')
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bo/phylo/<phylo_id>/feature/',
    view_func=phylo_feat,
    methods=['GET','POST']
)
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bo/phylo/<int:phylo_id>/feature/<int:cmt_id>',
    view_func=phylo,
    methods=['GET','PUT','DELETE']
)
