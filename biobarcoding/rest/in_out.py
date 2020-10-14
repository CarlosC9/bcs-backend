from flask import Blueprint

bp_io = Blueprint('io', __name__)

from flask import request, make_response, jsonify, send_file, abort
from flask.views import MethodView

import os
from biobarcoding.rest import bcs_api_base


class FileAPI(MethodView):
    """
    File Resource
    """
    def get(self, item):
        print(f'GET {request.path}\nDownloading file.')
        self._check_data()
        if item == 'sequence':
            from biobarcoding.services.sequences import export_sequences
            response, code = export_sequences(request.args['output'], sequence_id=self.sequence_id, organism_id=self.organism_id, analysis_id=self.analysis_id)
        # elif item == 'taxonomy':
        #     from biobarcoding.services.taxonomies import export_taxonomy
        #     response, code = export_taxonomy(request.args['output'], id=self.taxonomy_id)
        elif item == 'organism':
            from biobarcoding.services.organisms import export_organism
            response, code = export_organism(request.args['output'], self.organism_id)
        else:
            abort(make_response(jsonify({'status': 'failure', 'message': 'Method not available yet.'}), 405))
        if not code==200:
            abort(make_response(jsonify({'status': 'failure', 'message': response}), code))
        return send_file(response, mimetype='text/plain'), code


    def post(self, item):
        print(f'POST {request.path}\nUploading file.')
        self._check_data()
        msg = ''
        for key,file in request.files.items(multi=True):
            if file.filename == '':
                msg += 'status 409, message: File "" unknown. Not selected file.\n'
                continue
            file_cpy = self._make_file(file)
            resp = self._import_file(file_cpy, item)
            msg += f'FILE {os.path.basename(file_cpy)} RESULT: {resp}\n'
        responseObject = {
            'status': 'completed',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def _check_data(self):
        post_data = request.get_json()
        self.sequence_id = None
        self.taxonomy_id = None
        self.organism_id = None
        self.analysis_id = None
        if post_data:
            if 'sequence_id' in post_data:
                self.sequence_id = post_data['sequence_id']
            if 'taxonomy_id' in post_data:
                self.taxonomy_id = post_data['taxonomy_id']
            if 'organism_id' in post_data:
                self.organism_id = post_data['organism_id']
            if 'analysis_id' in post_data:
                self.analysis_id = post_data['analysis_id']
        print(f'JSON data: {post_data}')

    def _make_file(self, file):
        from werkzeug.utils import secure_filename
        file_path = os.path.join('/tmp', secure_filename(file.filename))
        file.save(file_path)
        return file_path

    def _import_file(self, filename, bioitem):
        if bioitem == 'sequence':
            from biobarcoding.services.sequences import import_sequences
            return import_sequences(filename, self.organism_id, self.analysis_id)
        elif bioitem == 'taxonomy':
            from biobarcoding.services.taxonomies import import_taxonomy
            return import_taxonomy(filename)
        elif bioitem == 'organism':
            abort(make_response(jsonify({'status': 'failure', 'message': 'Method not available yet.'}), 405))
        elif bioitem == 'ontology':
            abort(make_response(jsonify({'status': 'failure', 'message': 'Method not available yet.'}), 405))
        else:
            abort(make_response(jsonify({'status': 'failure', 'message': 'Method not available yet.'}), 405))

file = FileAPI.as_view('file_api')

bp_io.add_url_rule(
    bcs_api_base + '/io/<string:item>',
    view_func=file,
    methods=['GET','POST']
)
