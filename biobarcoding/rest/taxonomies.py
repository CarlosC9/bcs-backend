from flask import Blueprint

bp_taxonomies = Blueprint('bp_taxonomies', __name__)

from flask import request, make_response, jsonify, send_file
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


class TaxonomiesAPI(MethodView):
    """
    Taxonomies Resource
    """
    name = None
    comment = None

    def get(self, id=None):
        print(f'GET {request.path}\nGetting taxonomies {id}')
        self._check_data(request.args)
        if 'Accept' in request.headers and request.headers['Accept']=='text/ncbi':
            from biobarcoding.services.taxonomies import export_taxonomies
            response, code = export_taxonomies(id)
            return send_file(response, mimetype='text/ncbi'), code
        else:
            from biobarcoding.services.taxonomies import read_taxonomies
            response, code = read_taxonomies(id)
            return make_response(jsonify(response), code)


    def post(self):
        print(f'POST {request.path}\nCreating taxonomies')
        self._check_data(request.json)
        if 'Content-Type' in request.headers and request.headers['Content-Type']=='text/ncbi':
            response, code = self._import_files()
        else:
            from biobarcoding.services.taxonomies import create_taxonomies
            response, code = create_taxonomies(self.name, self.comment)
        return make_response(jsonify(response), code)


    def put(self, id):
        print(f'PUT {request.path}\nUpdating taxonomies {id}')
        self._check_data(request.json)
        from biobarcoding.services.taxonomies import update_taxonomies
        response, code = update_taxonomies(id, self.name, self.comment)
        return make_response(jsonify(response), code)


    def delete(self, id):
        print(f'DELETE {request.path}\nDeleting taxonomies {id}')
        from biobarcoding.services.taxonomies import delete_taxonomies
        response, code = delete_taxonomies(id)
        return make_response(jsonify(response), code)


    def _import_files(self):
        responses = []
        from biobarcoding.services.taxonomies import import_taxonomies
        for key,file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                response, code = import_taxonomies(file_cpy, self.name, self.comment)
                responses.append({'status':code,'message':response})
            except Exception as e:
                responses.append({'status':409,'message':'Could not import the file {file}.'})
        return responses, 207

    def _make_file(self, file):
        from werkzeug.utils import secure_filename
        file_path = os.path.join('/tmp', secure_filename(file.filename))
        file.save(file_path)
        return file_path


    def _check_data(self, data):
        if data:
            if 'name' in data:
                self.name = data['name']
            if 'comment' in data:
                self.comment = data['comment']
        print(f'DATA: {data}')


taxonomies_view = TaxonomiesAPI.as_view('api_taxonomies')
bp_taxonomies.add_url_rule(
    bcs_api_base + '/taxonomies/',
    view_func=taxonomies_view,
    methods=['GET','POST']
)
bp_taxonomies.add_url_rule(
    bcs_api_base + '/taxonomies/<int:id>',
    view_func=taxonomies_view,
    methods=['GET','PUT','DELETE']
)
