from flask import Blueprint

bp_alignments = Blueprint('msa', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


class AlignAPI(MethodView):
    """
    Align Resource
    """
    def get(self, id=None):
        msg = f'GET {request.path}\nGetting MSA {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def post(self):
        msg = f'POST {request.path}\nCreating MSA'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def put(self, id):
        msg = f'PUT {request.path}\nCreating MSA {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting MSA {id}'
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


class AlignFeatAPI(MethodView):
    """
    Alignment Feature Resource
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


msa = AlignAPI.as_view('msa_api')
bp_alignments.add_url_rule(
    bcs_api_base + '/bo/msa/',
    view_func=msa,
    methods=['GET','POST']
)
bp_alignments.add_url_rule(
    bcs_api_base + '/bo/msa/<int:msa_id>',
    view_func=msa,
    methods=['GET','PUT','DELETE']
)

msa_feat = AlignFeatAPI.as_view('msa_feat_api')
bp_alignments.add_url_rule(
    bcs_api_base + '/bo/msa/<msa_id>/feature/',
    view_func=msa_feat,
    methods=['GET','POST']
)
bp_alignments.add_url_rule(
    bcs_api_base + '/bo/msa/<int:msa_id>/feature/<int:cmt_id>',
    view_func=msa_feat,
    methods=['GET','PUT','DELETE']
)
