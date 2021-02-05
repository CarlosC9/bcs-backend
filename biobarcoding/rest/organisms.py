from flask import Blueprint, request, send_file
from flask.views import MethodView

bp_organisms = Blueprint('bp_organisms', __name__)

from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject


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

    @bcs_session(read_only=True)
    def get(self, id=None, format=None):
        print(f'GET {request.path}\nGetting organisms {id}')
        self._check_data(request.args)
        if format:
            from biobarcoding.services.organisms import export_organisms
            issues, content, status = export_organisms(id, format)
            return send_file(content, mimetype=f'text/{format}'), status
        from biobarcoding.services.organisms import read_organisms
        issues, content, status = read_organisms(id, self.genus, self.species, self.common_name, self.abbreviation, self.comment)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def post(self):
        print(f'POST {request.path}\nCreating organisms')
        self._check_data(request.json)
        from biobarcoding.services.organisms import create_organisms
        issues, content, status = create_organisms(self.genus, self.species, self.common_name, self.abbreviation, self.comment)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def put(self, id):
        print(f'PUT {request.path}\nUpdating organisms')
        self._check_data(request.json)
        from biobarcoding.services.organisms import update_organisms
        issues, content, status = update_organisms(id, self.genus, self.species, self.common_name, self.abbreviation, self.comment)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def delete(self, id=None):
        print(f'DELETE {request.path}\nDeleting organisms {id}')
        self._check_data(request.args)
        from biobarcoding.services.organisms import delete_organisms
        issues, content, status = delete_organisms(id, self.ids, self.genus, self.species, self.common_name, self.abbreviation, self.comment)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


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
bp_organisms.add_url_rule(
    bcs_api_base + '/organisms/<int:id>.<string:format>',
    view_func=organisms_view,
    methods=['GET','PUT','DELETE']
)
