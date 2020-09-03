from flask import Blueprint

bp_io = Blueprint('io', __name__)

from flask import request, make_response, jsonify, send_file
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


class FileAPI(MethodView):
    """
    File Resource
    """
    def get(self, item):
        print(f'GET {request.path}\nUploading file.')
        self._check_data()
        if item == 'sequence':
            from biobarcoding.services.sequences import export_sequences
            res = export_sequences(request.args['output'], self.organism_id)
        elif item == 'taxonomy':
            from biobarcoding.services.taxonomies import export_taxonomy
            res = export_taxonomy(request.files.items['file'])
        elif item == 'ontology':
            return make_response(jsonify({'status': 'failure', 'message': 'Method not available yet.'})), 405
        else:
            return make_response(jsonify({'status': 'failure', 'message': 'Method not available yet.'})), 405
        return send_file(res, mimetype='text/plain'), 200


    def post(self, item):
        print(f'POST {request.path}\nDownloading file.')
        self._check_data()

        msg = self._import(request, item)

        responseObject = {
            'status': 'completed',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def _import(self, request, item):
        msg = ''
        import os
        from werkzeug.utils import secure_filename
        for key,file in request.files.items(multi=True):
            if file.filename == '':
                msg += 'status 409, message: File "" unknown. Not selected file.\n'
                continue
            file_cpy = os.path.join('/tmp', secure_filename(file.filename))
            file.save(file_cpy)
            if item == 'sequence':
                from biobarcoding.services.sequences import import_sequences
                resp = import_sequences(file_cpy, self.organism_id)
            elif item == 'taxonomy':
                from biobarcoding.services.taxonomies import import_taxonomy
                resp = import_taxonomy(file_cpy)
            elif item == 'ontology':
                return make_response(jsonify({'status': 'failure', 'message': 'Method not available yet.'})), 405
            else:
                return make_response(jsonify({'status': 'failure', 'message': 'Method not available yet.'})), 405
            msg += f'FILE {file_cpy[5:]} RESULT: {resp}\n'
        return msg


    def _check_data(self):
        post_data = request.get_json()
        self.organism_id = None
        if 'organism_id' in post_data:
            self.organism_id = post_data['organism_id']
        print(f'JSON data: {post_data}')

file = FileAPI.as_view('file_api')

bp_io.add_url_rule(
    bcs_api_base + '/io/<string:item>',
    view_func=file,
    methods=['GET','POST']
)
