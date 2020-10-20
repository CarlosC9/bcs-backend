from flask import Blueprint

bp_ontologies = Blueprint('bp_ontologies', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


class OntologiesAPI(MethodView):
    """
    Ontologies Resource
    """
    name = None
    definition = None
    remote_url = None

    def get(self, id=None):
        print(f'GET {request.path}\nGetting ontologies {id}')
        # self._check_data(request.json)
        self._check_data(request.args)
        if 'Accept' in request.headers and request.headers['Accept']=='text/obo':
            from biobarcoding.services.ontologies import export_ontologies
            response, code = export_ontologies(id)
            return send_file(response, mimetype='text/obo'), code
        else:
            from biobarcoding.services.ontologies import read_ontologies
            response, code = read_ontologies(id, self.name)
            return make_response(jsonify(response), code)


    def post(self):
        print(f'POST {request.path}\nCreating ontologies')
        self._check_data(request.json)
        if 'Content-Type' in request.headers and request.headers['Content-Type']=='text/obo':
            response, code = self._import_files()
        else:
            from biobarcoding.services.ontologies import create_ontologies
            response, code = create_ontologies(self.name, definition = self.definition, remote_url = self.remote_url)
        return make_response(jsonify(response), code)


    def put(self, id):
        print(f'PUT {request.path}\nUpdating ontologies {id}')
        self._check_data(request.json)
        from biobarcoding.services.ontologies import update_ontologies
        response, code = update_ontologies(id,
            name=self.name,
            definition=self.definition,
            remote_url=self.remote_url,
            input_file=None)
        return make_response(jsonify(response), code)


    def delete(self, id):
        print(f'DELETE {request.path}\nDeleting ontologies {id}')
        # self._check_data(request.json)
        self._check_data(request.args)
        from biobarcoding.services.ontologies import delete_ontologies
        response, code = delete_ontologies(id)
        return make_response(jsonify(response), code)


    def _import_files(self):
        responses = []
        from biobarcoding.services.ontologies import import_ontologies
        for key,file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                response, code = import_ontologies(file_cpy, name = self.name, definition = self.definition)
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
            if 'name' in data and data['name']:
                self.name = data['name']
            if 'definition' in data and data['definition']:
                self.definition = data['definition']
            if 'remote_url' in data and data['remote_url']:
                self.remote_url = data['remote_url']
        print(f'DATA: {data}')


ontologies_view = OntologiesAPI.as_view('api_ontologies')
bp_ontologies.add_url_rule(
    bcs_api_base + '/ontologies/',
    view_func=ontologies_view,
    methods=['GET','POST']
)
bp_ontologies.add_url_rule(
    bcs_api_base + '/ontologies/<int:id>',
    view_func=ontologies_view,
    methods=['GET','PUT','DELETE']
)
