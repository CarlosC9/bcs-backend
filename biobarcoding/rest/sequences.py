from flask import Blueprint

from biobarcoding.authentication import bcs_session

bp_sequences = Blueprint('bp_sequences', __name__)

from flask import request, send_file
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base, ResponseObject, Issue, IType


class SequencesAPI(MethodView):
    """
    Sequences Resource
    """
    ids = None
    organism_id = None
    analysis_id = None
    phylotree_id = None
    format = 'fasta'

    @bcs_session(read_only=True)
    def get(self, id=None, format=None):
        print(f'GET {request.path}\nGetting sequences {id}')
        self._check_data(request.json)
        self._check_data(request.args)
        if format:
            from biobarcoding.services.sequences import export_sequences
            issues, content, status = export_sequences(id, self.ids, self.organism_id, self.analysis_id)
            return send_file(content, mimetype=f'text/{format}'), status
        from biobarcoding.services.sequences import read_sequences
        issues, content, status = read_sequences(id, self.ids, self.organism_id, self.analysis_id, self.phylotree_id)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def post(self):
        print(f'POST {request.path}\nCreating sequences')
        self._check_data(request.json)
        self._check_data(request.args)
        if request.files:
            issues, content, status = self._import_files()
        else:
            from biobarcoding.services.sequences import create_sequences
            issues, content, status = create_sequences(self.organism_id, self.analysis_id)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def put(self, id):
        print(f'PUT {request.path}\nCreating sequences {id}')
        self._check_data(request.json)
        from biobarcoding.services.sequences import update_sequences
        issues, content, status = update_sequences(id, self.organism, self.analysis)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def delete(self, id=None):
        print(f'DELETE {request.path}\nDeleting sequences {id}')
        self._check_data(request.json)
        self._check_data(request.args)
        from biobarcoding.services.sequences import delete_sequences
        issues, content, status = delete_sequences(id, self.ids, self.organism_id, self.analysis_id)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    def _import_files(self):
        from biobarcoding.services.sequences import import_sequences
        issues, content = [], []
        for key,file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                i, c, s = import_sequences(file_cpy, self.organism_id, self.analysis_id, self.format)
            except Exception as e:
                print(e)
                i, c = [Issue(IType.ERROR, f'Could not import the file {file}.')], {}
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
            if 'organism_id' in data and data['organism_id']:
                self.organism_id = data['organism_id']
            if 'analysis_id' in data and data['analysis_id']:
                self.analysis_id = data['analysis_id']
            if 'phylotree_id' in data and data['phylotree_id']:
                self.phylotree_id = data['phylotree_id']
            if 'format' in data and data['format']:
                self.format = data['format']
        print(f'DATA: {data}')


sequences_view = SequencesAPI.as_view('api_sequences')
bp_sequences.add_url_rule(
    bcs_api_base + '/bos/sequences/',
    view_func=sequences_view,
    methods=['GET','POST','DELETE']
)
bp_sequences.add_url_rule(
    bcs_api_base + '/bos/sequences/<string:id>',
    view_func=sequences_view,
    methods=['GET','PUT','DELETE']
)
bp_sequences.add_url_rule(
    bcs_api_base + '/bos/sequences.<string:format>',
    view_func=sequences_view,
    methods=['GET']
)
bp_sequences.add_url_rule(
    bcs_api_base + '/bos/sequences/<string:id>.<string:format>',
    view_func=sequences_view,
    methods=['GET']
)


class SeqFeatAPI(MethodView):
    """
    Sequences Feature Resource
    """
    def get(self, seq_id, cmt_id=None):
        print(f'GET {request.path}\nGetting comments {seq_id} {cmt_id}')
        return ResponseObject(content={'status':'success','message':'READ: sequence comments, dummy complete'}, status=200).get_response()


    def post(self, seq_id):
        print(f'POST {request.path}\nCreating comments')
        return ResponseObject(content={'status':'success','message':'CREATE: sequence comments, dummy complete'}, status=200).get_response()


    def put(self, seq_id, cmt_id):
        print(f'PUT {request.path}\nCreating comments {seq_id} {cmt_id}')
        return ResponseObject(content={'status':'success','message':'UPDATE: sequence comments, dummy complete'}, status=200).get_response()


    def delete(self, seq_id, cmt_id=None):
        print(f'DELETE {request.path}\nDeleting comments {seq_id} {cmt_id}')
        return ResponseObject(content={'status':'success','message':'DELETE: sequence comments, dummy complete'}, status=200).get_response()


    def _check_data(self, data):
        if data:
            if 'id' in data and data['id']:
                self.ids = data['id']
            if 'comment' in data and data['comment']:
                self.comment = data['comment']
            if 'range' in data and data['range']:
                self.range = data['range']
        print(f'DATA: {data}')


seq_feat_view = SeqFeatAPI.as_view('api_seq_feat')
bp_sequences.add_url_rule(
    bcs_api_base + '/bos/sequences/<int:seq_id>/features/',
    view_func=seq_feat_view,
    methods=['GET','POST']
)
bp_sequences.add_url_rule(
    bcs_api_base + '/bos/sequences/<int:seq_id>/features/<int:cmt_id>',
    view_func=seq_feat_view,
    methods=['GET','PUT','DELETE']
)
