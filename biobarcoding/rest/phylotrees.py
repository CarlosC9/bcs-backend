from flask import Blueprint, request, send_file
from flask.views import MethodView

bp_phylotrees = Blueprint('bp_phylotrees', __name__)

from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject, Issue, IType


class PhyloAPI(MethodView):
    """
    Phylo Resource
    """
    ids = None
    analysis_id = None
    name = None
    comment = None
    feature_id = None
    format = None

    @bcs_session(read_only=True)
    def get(self, id=None, format=None):
        print(f'GET {request.path}\nGetting phylotrees {id}')
        self._check_data(request.json)
        self._check_data(request.args)
        if format:
            from biobarcoding.services.phylotrees import export_phylotrees
            issues, content, status = export_phylotrees(id, format)
            return send_file(content, mimetype=f'text/{format}'), status
        from biobarcoding.services.phylotrees import read_phylotrees
        issues, content, status = read_phylotrees(id, self.analysis_id, self.name, self.comment, self.feature_id)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def post(self):
        print(f'POST {request.path}\nCreating phylotrees')
        self._check_data(request.json)
        self._check_data(request.args)
        if request.files:
            issues, content, status = self._import_files()
        else:
            from biobarcoding.services.phylotrees import create_phylotrees
            issues, content, status = create_phylotrees(self.name, self.comment)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def put(self, id):
        print(f'PUT {request.path}\nUpdating phylotrees {id}')
        self._check_data(request.json)
        from biobarcoding.services.phylotrees import update_phylotrees
        issues, content, status = update_phylotrees(id, self.name, self.comment)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def delete(self, id=None):
        print(f'DELETE {request.path}\nDeleting phylotrees {id}')
        self._check_data(request.args)
        from biobarcoding.services.phylotrees import delete_phylotrees
        issues, content, status = delete_phylotrees(id, self.ids)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    def _import_files(self):
        issues, content = [], []
        self.format = self.format or 'newick'
        from biobarcoding.services.phylotrees import import_phylotrees
        for key,file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                i, c, s = import_phylotrees(file_cpy, analysis_id=self.analysis_id, name=self.name, comment=self.comment)
            except Exception as e:
                print(e)
                i, c = [Issue(IType.ERROR, f'Could not import the file {file}.')], None
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
            if 'analysis_id' in data:
                self.analysis_id = data['analysis_id']
            if 'name' in data:
                self.name = data['name']
            if 'comment' in data:
                self.comment = data['comment']
            if 'feature_id' in data:
                self.feature_id = data['feature_id']
            if 'format' in data and data['format']:
                self.format = data['format']
        print(f'DATA: {data}')


phylotrees_view = PhyloAPI.as_view('api_phylotrees')
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bos/phylotrees/',
    view_func=phylotrees_view,
    methods=['GET','POST','DELETE']
)
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bos/phylotrees/<int:id>',
    view_func=phylotrees_view,
    methods=['GET','PUT','DELETE']
)
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bos/phylotrees.<string:format>',
    view_func=phylotrees_view,
    methods=['GET','DELETE']
)
bp_phylotrees.add_url_rule(
    bcs_api_base + '/bos/phylotrees/<int:id>.<string:format>',
    view_func=phylotrees_view,
    methods=['GET','POST','PUT','DELETE']
)


class PhyloFeatAPI(MethodView):
    """
    Phylogenetic Tree Feature Resource
    """
    @bcs_session(read_only=True)
    def get(self, phylo_id, cmt_id=None):
        print(f'GET {request.path}\nGetting comment {phylo_id}:{cmt_id}')
        issues = [Issue(IType.WARNING, f'GET phylogeny comment: dummy completed.')]
        return ResponseObject(issues=issues).get_response()


    @bcs_session()
    def post(self, phylo_id):
        print(f'POST {request.path}\nCreating comment {phylo_id}')
        issues = [Issue(IType.WARNING, f'POST phylogeny comment: dummy completed.')]
        return ResponseObject(issues=issues).get_response()


    @bcs_session()
    def put(self, phylo_id, cmt_id):
        print(f'PUT {request.path}\nUpdating comment {phylo_id} {cmt_id}')
        issues = [Issue(IType.WARNING, f'PUT phylogeny comment: dummy completed.')]
        return ResponseObject(issues=issues).get_response()


    @bcs_session()
    def delete(self, phylo_id, cmt_id=None):
        print(f'DELETE {request.path}\nDeleting comment {phylo_id} {cmt_id}')
        issues = [Issue(IType.WARNING, f'DELETE phylogeny comment: dummy completed.')]
        return ResponseObject(issues=issues).get_response()


    def _check_data(self):
        post_data = request.json
        print(f'JSON data: {post_data}')


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
