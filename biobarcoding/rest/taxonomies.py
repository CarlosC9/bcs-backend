from flask import Blueprint, request, send_file
from flask.views import MethodView

bp_taxonomies = Blueprint('bp_taxonomies', __name__)

from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject, Issue, IType


class TaxonomiesAPI(MethodView):
    """
    Taxonomies Resource
    """
    ids = None
    name = None
    comment = None

    @bcs_session(read_only=True)
    def get(self, id=None, format=None):
        print(f'GET {request.path}\nGetting taxonomies {id}')
        self._check_data(request.args)
        if format:
            from biobarcoding.services.taxonomies import export_taxonomies
            issues, content, status = export_taxonomies(id, self.ids, format)
            return send_file(content, mimetype='text/{format}'), status
        from biobarcoding.services.taxonomies import read_taxonomies
        issues, content, status = read_taxonomies(id, self.name, self.comment)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def post(self):
        print(f'POST {request.path}\nCreating taxonomies')
        self._check_data(request.args)
        self._check_data(request.json)
        if request.files:
            issues, content, status = self._import_files()
        else:
            from biobarcoding.services.taxonomies import create_taxonomies
            issues, content, status = create_taxonomies(self.name, self.comment)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def put(self, id):
        print(f'PUT {request.path}\nUpdating taxonomies {id}')
        self._check_data(request.json)
        from biobarcoding.services.taxonomies import update_taxonomies
        issues, content, status = update_taxonomies(id, self.name, self.comment)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def delete(self, id=None):
        print(f'DELETE {request.path}\nDeleting taxonomies {id}')
        self._check_data(request.args)
        from biobarcoding.services.taxonomies import delete_taxonomies
        issues, content, status = delete_taxonomies(id, self.ids)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    def _import_files(self):
        issues, content = [], []
        from biobarcoding.services.taxonomies import import_taxonomies
        for key,file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                i, c, s = import_taxonomies(file_cpy, self.name, self.comment)
            except Exception as e:
                print(e)
                i, c = Issue(IType.ERROR, f'Could not import the file {file}.'), None
            issues+=i
            content.append(c)
        return issues, content, 207


    def _make_file(self, file):
        import os
        from werkzeug.utils import secure_filename
        file_path = os.path.join('/tmp', secure_filename(file.filename))
        file.save(file_path)
        return file_path


    def _check_data(self, data):
        if data:
            if 'id' in data and data['id']:
                self.ids = data.getlist('id')
            if 'name' in data and data['name']:
                self.name = data['name']
            if 'comment' in data and data['comment']:
                self.comment = data['comment']
        print(f'DATA: {data}')


taxonomies_view = TaxonomiesAPI.as_view('api_taxonomies')
bp_taxonomies.add_url_rule(
    bcs_api_base + '/taxonomies/',
    view_func=taxonomies_view,
    methods=['GET','POST','DELETE']
)
bp_taxonomies.add_url_rule(
    bcs_api_base + '/taxonomies/<int:id>',
    view_func=taxonomies_view,
    methods=['GET','PUT','DELETE']
)
bp_taxonomies.add_url_rule(
    bcs_api_base + '/taxonomies/<int:id>.<string:format>',
    view_func=taxonomies_view,
    methods=['GET','PUT','DELETE']
)
