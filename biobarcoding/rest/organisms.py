from flask import Blueprint

bp_organisms = Blueprint('bp_organisms', __name__)

from flask import request, make_response, jsonify, send_file
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


# Organisms API

class OrganismsAPI(MethodView):
    """
    Organisms Resource
    """
    ids = None
    genus = None
    species = None
    common_name = None
    abbreviation = None
    comment = None

    def get(self, id=None):
        print(f'GET {request.path}\nGetting organisms {id}')
        self._check_data(request.args)
        if 'Accept' in request.headers and request.headers['Accept']=='text/genbank':
            from biobarcoding.services.organisms import export_organisms
            response, code = export_organisms(id)
            return send_file(response, mimetype='text/genbank'), code
        else:
            from biobarcoding.services.organisms import read_organisms
            response, code = read_organisms(id, self.genus, self.species, self.common_name, self.abbreviation, self.comment)
            return make_response(jsonify(response), code)


    def post(self):
        print(f'POST {request.path}\nCreating organisms')
        self._check_data(request.json)
        from biobarcoding.services.organisms import create_organisms
        response, code = create_organisms(self.genus, self.species, self.common_name, self.abbreviation, self.comment)
        return make_response(jsonify(response), code)


    def put(self, id):
        print(f'PUT {request.path}\nUpdating organisms')
        self._check_data(request.json)
        from biobarcoding.services.organisms import update_organisms
        response, code = update_organisms(id, self.genus, self.species, self.common_name, self.abbreviation, self.comment)
        return make_response(jsonify(response), code)


    def delete(self, id=None):
        print(f'DELETE {request.path}\nDeleting organisms {id}')
        self._check_data(request.args)
        from biobarcoding.services.organisms import delete_organisms
        response, code = delete_organisms(id, self.ids, self.genus, self.species, self.common_name, self.abbreviation, self.comment)
        return make_response(jsonify(response), code)


    def _check_data(self, data):
        if data:
            if 'id' in data and data['id']:
                self.ids = data.getlist('id')
            if 'genus' in data and data['genus']:
                self.genus = data['genus']
            if 'species' in data and data['species']:
                self.species = data['species']
            if 'common_name' in data and data['common_name']:
                self.common_name = data['common_name']
            if 'abbreviation' in data and data['abbreviation']:
                self.abbreviation = data['abbreviation']
            if 'comment' in data and data['comment']:
                self.comment = data['comment']
        print(f'DATA: {data}')


organisms_view = OrganismsAPI.as_view('api_organisms')
bp_organisms.add_url_rule(
    bcs_api_base + '/organisms/',
    view_func=organisms_view,
    methods=['GET','POST','DELETE']
)
bp_organisms.add_url_rule(
    bcs_api_base + '/organisms/<int:id>',
    view_func=organisms_view,
    methods=['GET','PUT','DELETE']
)
