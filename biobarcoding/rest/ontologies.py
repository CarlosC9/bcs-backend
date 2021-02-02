from flask import Blueprint, request, send_file
from flask.views import MethodView

bp_ontologies = Blueprint('bp_ontologies', __name__)

from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject, Issue, IType


class OntologiesAPI(MethodView):
    """
    Ontologies Resource
    """
    name = None
    definition = None
    remote_url = None
    format = None

    @bcs_session(read_only=True)
    def get(self, id=None, format=None):
        print(f'GET {request.path}\nGetting ontologies {id}')
        self._check_data(request.json)
        self._check_data(dict(request.values))
        if format or self.format:
            format = format or self.format
            from biobarcoding.services.ontologies import export_ontologies
            issues, content, status = export_ontologies(id, self.name, self.definition, self.remote_url, format)
            return send_file(content, mimetype=f'text/{format}'), status
        from biobarcoding.services.ontologies import read_ontologies
        issues, content, status = read_ontologies(id, self.name, self.definition, self.remote_url)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def post(self):
        print(f'POST {request.path}\nCreating ontologies')
        if request.files:
            issues, content, status = self._import_files()
        else:
            self._check_data(request.json)
            self._check_data(dict(request.values))
            from biobarcoding.services.ontologies import create_ontologies
            issues, content, status = create_ontologies(self.name, definition=self.definition, remote_url=self.remote_url)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def put(self, id):
        print(f'PUT {request.path}\nUpdating ontologies {id}')
        self._check_data(request.json)
        self._check_data(dict(request.values))
        from biobarcoding.services.ontologies import update_ontologies
        issues, content, status = update_ontologies(id,
                                                    name=self.name,
                                                    definition=self.definition,
                                                    remote_url=self.remote_url,
                                                    input_file=None)
        return ResponseObject(content=content, issues=issues, status=status).get_response()



    @bcs_session()
    def delete(self, id):
        print(f'DELETE {request.path}\nDeleting ontologies {id}')
        self._check_data(request.json)
        self._check_data(dict(request.values))
        from biobarcoding.services.ontologies import delete_ontologies
        issues, content, status = delete_ontologies(id)
        return ResponseObject(content=content, issues=issues, status=status).get_response()



    def _import_files(self):
        issues, content = [], []
        self.format = self.format or 'obo'
        from biobarcoding.services.ontologies import import_ontologies
        for key,file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                i, c, s = import_ontologies(file_cpy, format=self.format, name=self.name, definition=self.definition, remote_url=self.remote_url)
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
bp_ontologies.add_url_rule(
    bcs_api_base + '/ontologies/<int:id>.<string:format>',
    view_func=ontologies_view,
    methods=['GET','PUT','DELETE']
)


class CvtermsAPI(MethodView):
    """
    Cvterms Resource
    """
    feature_id = None
    analysis_id = None
    phylotree_id = None

    @bcs_session(read_only=True)
    def get(self, cv_id=None, cvterm_id=None):
        print(f'GET {request.path}\nGetting ontology terms {id}')
        self._check_data(request.json)
        self._check_data(dict(request.values))
        from biobarcoding.services.ontologies import read_cvterms
        issues, content, status = read_cvterms(cv_id, cvterm_id, self.feature_id, self.analysis_id, self.phylotree_id)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    def _check_data(self, data):
        if data:
            if 'feature_id' in data and data['feature_id']:
                self.feature_id = data['feature_id']
            if 'analysis_id' in data and data['analysis_id']:
                self.analysis_id = data['analysis_id']
            if 'phylotree_id' in data and data['phylotree_id']:
                self.phylotree_id = data['phylotree_id']
        print(f'DATA: {data}')

cvterms_view = CvtermsAPI.as_view('api_cvterms')
bp_ontologies.add_url_rule(
    bcs_api_base + '/ontologies/terms/',
    view_func=cvterms_view,
    methods=['GET']
)
bp_ontologies.add_url_rule(
    bcs_api_base + '/ontologies/<int:cv_id>/terms/',
    view_func=cvterms_view,
    methods=['GET']
)
bp_ontologies.add_url_rule(
    bcs_api_base + '/ontologies/terms/<int:cvterm_id>',
    view_func=cvterms_view,
    methods=['GET']
)
bp_ontologies.add_url_rule(
    bcs_api_base + '/ontologies/<int:cv_id>/terms/<int:cvterm_id>',
    view_func=cvterms_view,
    methods=['GET']
)
