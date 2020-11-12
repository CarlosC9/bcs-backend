from flask import Blueprint

bp_alignments = Blueprint('bp_alignments', __name__)

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

        responseObject = []
        return make_response(jsonify(responseObject), 200)


    def post(self):
        msg = f'POST {request.path}\nCreating MSA'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject), 200)


    def put(self, id):
        msg = f'PUT {request.path}\nCreating MSA {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject), 200)


    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting MSA {id}'
        print(msg)
        self._check_data()

        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject), 200)


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
        return make_response(jsonify(responseObject), 200)


    def post(self):
        msg = f'POST {request.path}\nCreating comment'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject), 200)


    def put(self, id):
        msg = f'PUT {request.path}\nCreating comment {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject), 200)


    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting comment {id}'
        print(msg)
        self._check_data()

        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject), 200)


    def _check_data(self):

        post_data = request.get_json()
        print(f'JSON data: {post_data}')


alignment_view = AlignAPI.as_view('api_alignments')
bp_alignments.add_url_rule(
    bcs_api_base + '/bos/alignments/',
    view_func=alignment_view,
    methods=['GET','POST']
)
bp_alignments.add_url_rule(
    bcs_api_base + '/bos/alignments/<int:alignment_id>',
    view_func=alignment_view,
    methods=['GET','PUT','DELETE']
)

alignment_feat_view = AlignFeatAPI.as_view('api_alignment_feat')
bp_alignments.add_url_rule(
    bcs_api_base + '/bos/alignments/<alignment_id>/features/',
    view_func=alignment_feat_view,
    methods=['GET','POST']
)
bp_alignments.add_url_rule(
    bcs_api_base + '/bos/alignments/<int:alignment_id>/features/<int:cmt_id>',
    view_func=alignment_feat_view,
    methods=['GET','PUT','DELETE']
)
