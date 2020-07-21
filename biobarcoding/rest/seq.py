from flask import Blueprint

bp_seq = Blueprint('seq', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

class SeqAPI(MethodView):
    """
    Sequence Resource
    """
    def get(self, id=None):
        msg = f'GET {request.path}\nGetting sequence {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def post(self):
        msg = f'POST {request.path}\nCreating sequence'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def put(self, id):
        msg = f'PUT {request.path}\nCreating sequence {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting sequence {id}'
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


seq = SeqAPI.as_view('seq_api')
bp_seq.add_url_rule(
    '/bo/sequence',
    view_func=seq,
    methods=['GET','POST']
)
bp_seq.add_url_rule(
    '/bo/sequence/<int:id>',
    view_func=seq,
    methods=['GET','PUT','DELETE']
)


class SeqFeatAPI(MethodView):
    """
    Sequence Feature Resource
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


seq_feat = SeqFeatAPI.as_view('seq_feat_api')
bp_seq.add_url_rule(
    '/bo/sequence/<seq_id>/feature/',
    view_func=seq_feat,
    methods=['GET','POST']
)
bp_seq.add_url_rule(
    '/bo/sequence/<int:seq_id>/feature/<int:cmt_id>',
    view_func=seq_feat,
    methods=['GET','PUT','DELETE']
)
