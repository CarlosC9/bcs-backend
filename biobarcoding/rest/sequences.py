from flask import Blueprint

bp_sequences = Blueprint('sequences', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


class SequencesAPI(MethodView):
    """
    Sequences Resource
    """
    ids = None
    organism_id = None
    analysis_id = None

    def get(self, id=None):
        print(f'GET {request.path}\nGetting sequences {id}')
        self._check_data(request.get_json())
        self._check_data(request.args.to_dict())
        if 'Accept' in request.headers and request.headers['Accept']=='text/fasta':
            from biobarcoding.services.sequences import export_sequences
            response, code = export_sequences(id, self.organism_id, self.analysis_id)
            return send_file(response, mimetype='text/fasta'), code
        else:
            from biobarcoding.services.sequences import read_sequences
            response, code = read_sequences(id, self.organism_id, self.analysis_id)
            return make_response(response, code)


    def post(self):
        print(f'POST {request.path}\nCreating sequences')
        self._check_data(request.get_json())
        if 'Content-Type' in request.headers and request.headers['Content-Type']=='text/fasta':
            response, code = self._import_files()
        else:
            from biobarcoding.services.sequences import create_sequences
            response, code = create_sequences(self.organism_id, self.analysis_id)
        return make_response(response, code)


    def put(self, id):
        print(f'PUT {request.path}\nCreating sequences {id}')
        self._check_data(request.get_json())
        from biobarcoding.services.sequences import update_sequences
        response, code = update_sequences(id, self.organism, self.analysis)
        return make_response(response, code)


    def delete(self, id=None):
        print(f'DELETE {request.path}\nDeleting sequences {id}')
        self._check_data(request.get_json())
        self._check_data(request.args.to_dict())
        from biobarcoding.services.sequences import delete_sequences
        response, code = delete_sequences(id, self.organism_id, self.analysis_id)
        return make_response(response, code)


    def _import_files(self):
        responses = []
        from biobarcoding.services.sequences import import_sequences
        for key,file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                response, code = import_sequences(file_cpy, self.organism_id, self.analysis_id)
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
            if 'id' in data and data['id']:
                self.ids = data['id']
            if 'organism_id' in data and data['organism_id']:
                self.organism_id = data['organism_id']
            if 'analysis_id' in data and data['analysis_id']:
                self.analysis_id = data['analysis_id']
        print(f'DATA: {data}')


sequences = SequencesAPI.as_view('sequences_api')
bp_sequences.add_url_rule(
    bcs_api_base + '/bos/sequences/',
    view_func=sequences,
    methods=['GET','POST','DELETE']
)
bp_sequences.add_url_rule(
    bcs_api_base + '/bos/sequences/<int:id>',
    view_func=sequences,
    methods=['GET','PUT','DELETE']
)


class SeqFeatAPI(MethodView):
    """
    Sequences Feature Resource
    """
    def get(self, seq_id, cmt_id=None):
        print(f'GET {request.path}\nGetting comments {seq_id} {cmt_id}')
        # self._check_data(request.args.to_dict())
        return make_response({'status':'success','message':'dummy complete'}, 200)


    def post(self, seq_id):
        print(f'POST {request.path}\nCreating comments')
        # self._check_data(request.args.to_dict())
        return make_response({'status':'success','message':'dummy complete'}, 200)


    def put(self, seq_id, cmt_id):
        print(f'PUT {request.path}\nCreating comments {seq_id} {cmt_id}')
        # self._check_data(request.args.to_dict())
        return make_response({'status':'success','message':'dummy complete'}, 200)


    def delete(self, seq_id, cmt_id=None):
        print(f'DELETE {request.path}\nDeleting comments {seq_id} {cmt_id}')
        # self._check_data(request.args.to_dict())
        return make_response({'status':'success','message':'dummy complete'}, 200)


    def _check_data(self, data):
        if data:
            if 'id' in data and data['id']:
                self.ids = data['id']
            if 'comment' in data and data['comment']:
                self.comment = data['comment']
            if 'range' in data and data['range']:
                self.range = data['range']
        print(f'DATA: {data}')


seq_feat = SeqFeatAPI.as_view('seq_feat_api')
bp_sequences.add_url_rule(
    bcs_api_base + '/bos/sequences/<int:seq_id>/features/',
    view_func=seq_feat,
    methods=['GET','POST']
)
bp_sequences.add_url_rule(
    bcs_api_base + '/bos/sequences/<int:seq_id>/features/<int:cmt_id>',
    view_func=seq_feat,
    methods=['GET','PUT','DELETE']
)
