from flask import Blueprint

bp_onto = Blueprint('ontology', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

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
bp_onto.add_url_rule(
    '/bo/ontology',
    view_func=onto,
    methods=['GET','POST']
)
bp_onto.add_url_rule(
    '/bo/ontology/<int:id>',
    view_func=onto,
    methods=['GET','DELETE']
)
