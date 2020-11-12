from flask import Blueprint

bp_phylotrees = Blueprint('bp_phylotrees', __name__)

from flask import request, make_response, jsonify, send_file
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


class PhyloAPI(MethodView):
    """
    Phylo Resource
    """
    analysis_id = None
    name = None
    comment = None

    def get(self, id=None):
        print(f'GET {request.path}\nGetting phylotrees {id}')
        self._check_data(request.json)
        if 'Accept' in request.headers and request.headers['Accept']=='text/newick':
            from biobarcoding.services.phylotrees import export_phylotrees
            response, code = export_phylotrees(id)
            return send_file(response, mimetype='text/newick'), code
        else:
            from biobarcoding.services.phylotrees import read_phylotrees
            response, code = read_phylotrees(id)
        return make_response(jsonify(response), code)


    def post(self):
        print(f'POST {request.path}\nCreating phylotrees')
        self._check_data(request.json)
        self._check_data(request.args)
        if request.files:
            response, code = self._import_files()
        else:
            from biobarcoding.services.phylotrees import create_phylotrees
            response, code = create_phylotrees(self.name, self.comment)
        return make_response(jsonify(response), code)


    def put(self, id):
        print(f'PUT {request.path}\nUpdating phylotrees {id}')
        self._check_data(request.json)
        from biobarcoding.services.phylotrees import update_phylotrees
        response, code = update_phylotrees(id, self.name, self.comment)
        return make_response(jsonify(response), code)


    def delete(self, id):
        print(f'DELETE {request.path}\nDeleting phylotrees {id}')
        from biobarcoding.services.phylotrees import delete_phylotrees
        response, code = delete_phylotrees(id)
        return make_response(jsonify(response), code)


    def _import_files(self):
        responses = []
        from biobarcoding.services.phylotrees import import_phylotrees
        for key,file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                response, code = import_phylotrees(file_cpy, analysis_id=self.analysis_id, name=self.name, comment=self.comment)
                responses.append(response)
            except Exception as e:
                print(e)
                responses.append({'status':409,'message':f'Could not import the file {file}.'})
        return responses, 207

    def _make_file(self, file):
        import os
        from werkzeug.utils import secure_filename
        file_path = os.path.join('/tmp', secure_filename(file.filename))
        file.save(file_path)
        return file_path


    def _check_data(self, data):
        if data:
            if 'analysis_id' in data:
                self.analysis_id = data['analysis_id']
            if 'name' in data:
                self.name = data['name']
            if 'comment' in data:
                self.comment = data['comment']
        print(f'DATA: {data}')


class PhyloFeatAPI(MethodView):
    """
    Phylogenetic Tree Feature Resource
    """
    def get(self, id=None):
        msg = f'GET {request.path}\nGetting comment {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def post(self):
        msg = f'POST {request.path}\nCreating comment'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def put(self, id):
        msg = f'PUT {request.path}\nCreating comment {id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting comment {id}'
        print(msg)
        self._check_data()

        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def _check_data(self):

        post_data = request.json
        print(f'JSON data: {post_data}')


phylotrees_view = PhyloAPI.as_view('api_phylotrees')
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bos/phylotrees/',
    view_func=phylotrees_view,
    methods=['GET','POST']
)
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bos/phylotrees/<int:phylo_id>',
    view_func=phylotrees_view,
    methods=['GET','PUT','DELETE']
)

phylo_feat_view = PhyloFeatAPI.as_view('api_phylo_feat')
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bos/phylotrees/<phylo_id>/features/',
    view_func=phylo_feat_view,
    methods=['GET','POST']
)
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bos/phylotrees/<int:phylo_id>/features/<int:cmt_id>',
    view_func=phylo_feat_view,
    methods=['GET','PUT','DELETE']
)
