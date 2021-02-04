from flask import Blueprint, request, send_file
from flask.views import MethodView

bp_alignments = Blueprint('bp_alignments', __name__)

from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject, Issue, IType


class AlignAPI(MethodView):
    """
    Align Resource
    """
    analysis_id = None
    ids = None
    format = None
    program = None
    programversion = None
    name = None
    sourcename = None
    description = None
    algorithm = None
    sourceversion = None
    sourceuri = None
    timeexecuted = None
    feature_id = None
    format = None

    @bcs_session(read_only=True)
    def get(self, alignment_id=None, format=None):
        print(f'GET {request.path}\nGetting alignments {alignment_id}')
        self._check_data(request.args)
        if format or self.format:
            format = format or self.format
            from biobarcoding.services.alignments import export_alignments
            issues, content, status = export_alignments(analysis_id=alignment_id, format=format)
            return send_file(content, mimetype=f'text/{format}'), status
        from biobarcoding.services.alignments import read_alignments
        issues, content, status = read_alignments(
            alignment_id=alignment_id,
            ids=self.ids,
            name=self.name,
            program=self.program,
            programversion=self.programversion,
            algorithm=self.algorithm,
            sourcename=self.sourcename,
            sourceversion=self.sourceversion,
            sourceuri=self.sourceuri,
            description=self.description,
            feature_id=self.feature_id)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def post(self):
        print(f'POST {request.path}\nCreating alignments')
        self._check_data(request.json)
        self._check_data(request.args)
        if request.files:
            issues, content, status = self._import_files()
        else:
            from biobarcoding.services.alignments import create_alignments
            issues, content, status = create_alignments(
                program=self.program or 'Multiple Sequence Alignment',
                programversion=self.programversion or 'unknown',
                name=self.name,
                sourcename=self.sourcename,
                description=self.description,
                algorithm=self.algorithm,
                sourceversion=self.sourceversion,
                sourceuri=self.sourceuri,
                timeexecuted=self.timeexecuted)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def put(self, alignment_id):
        print(f'PUT {request.path}\nUpdating alignments')
        self._check_data(request.json)
        from biobarcoding.services.alignments import update_alignments
        issues, content, status = update_alignments(alignment_id,
                                           program=self.program,
                                           programversion=self.programversion,
                                           name=self.name,
                                           sourcename=self.sourcename,
                                           description=self.description,
                                           algorithm=self.algorithm,
                                           sourceversion=self.sourceversion,
                                           sourceuri=self.sourceuri,
                                           timeexecuted=self.timeexecuted)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def delete(self, alignment_id=None):
        print(f'DELETE {request.path}\nDeleting alignments {alignment_id}')
        self._check_data(request.args)
        from biobarcoding.services.alignments import delete_alignments
        issues, content, status = delete_alignments(alignment_id,
                                           ids=self.ids,
                                           name=self.name,
                                           program=self.program,
                                           programversion=self.programversion,
                                           algorithm=self.algorithm,
                                           sourcename=self.sourcename,
                                           sourceversion=self.sourceversion,
                                           sourceuri=self.sourceuri,
                                           description=self.description)
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    def _import_files(self):
        issues, content = [], []
        self.format = self.format or 'fasta'
        from biobarcoding.services.alignments import import_alignments
        for key, file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                i, c, s = import_alignments(file_cpy, format=self.format,
                                                   analysis_id=self.analysis_id,
                                                   program=self.program or 'Multiple Sequence Alignment',
                                                   programversion=self.programversion or f'(Imported {self.format})',
                                                   name=self.name,
                                                   sourcename=self.sourcename,
                                                   description=self.description,
                                                   algorithm=self.algorithm,
                                                   sourceversion=self.sourceversion,
                                                   sourceuri=self.sourceuri,
                                                   timeexecuted=self.timeexecuted)
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
            if 'analysis_id' in data and data['analysis_id']:
                self.analysis_id = data['analysis_id']
            if 'id' in data and data['id']:
                self.ids = data.getlist('id')
            if 'program' in data and data['program']:
                self.program = data['program']
            if 'format' in data and data['format']:
                self.format = data['format']
                programversion = f'(Imported {format})'
            if 'programversion' in data and data['programversion']:
                self.programversion = data['programversion']
            if 'name' in data and data['name']:
                self.name = data['name']
            if 'sourcename' in data and data['sourcename']:
                self.sourcename = data['sourcename']
            if 'description' in data and data['description']:
                self.description = data['description']
            if 'algorithm' in data and data['algorithm']:
                self.algorithm = data['algorithm']
            if 'sourceversion' in data and data['sourceversion']:
                self.sourceversion = data['sourceversion']
            if 'sourceuri' in data and data['sourceuri']:
                self.sourceuri = data['sourceuri']
            if 'timeexecuted' in data and data['timeexecuted']:
                self.timeexecuted = data['timeexecuted']
            if 'feature_id' in data and data['feature_id']:
                self.feature_id = data['feature_id']
        print(f'DATA: {data}')


alignment_view = AlignAPI.as_view('api_alignments')
bp_alignments.add_url_rule(
    bcs_api_base + '/bos/alignments/',
    view_func=alignment_view,
    methods=['GET', 'POST', 'DELETE']
)
bp_alignments.add_url_rule(
    bcs_api_base + '/bos/alignments/<int:alignment_id>',
    view_func=alignment_view,
    methods=['GET', 'PUT', 'DELETE']
)
bp_alignments.add_url_rule(
    bcs_api_base + '/bos/alignments/<int:alignment_id>.<string:format>',
    view_func=alignment_view,
    methods=['GET', 'PUT', 'DELETE']
)

class AlignFeatAPI(MethodView):
    """
    Alignment Feature Resource
    """

    def get(self, cmt_id=None):
        print(f'GET {request.path}\nGetting comment {cmt_id}')
        self._check_data()
        return 'dummy response', 200

    def post(self):
        print(f'POST {request.path}\nCreating comment')
        self._check_data()
        return 'dummy response', 200

    def put(self, cmt_id):
        print(f'PUT {request.path}\nCreating comment {cmt_id}')
        self._check_data()
        return 'dummy response', 200

    def delete(self, cmt_id):
        print(f'DELETE {request.path}\nDeleting comment {cmt_id}')
        self._check_data()
        return 'dummy response', 200

    def _check_data(self):
        post_data = request.get_json()
        print(f'JSON data: {post_data}')


alignment_feat_view = AlignFeatAPI.as_view('api_alignment_feat')
bp_alignments.add_url_rule(
    bcs_api_base + '/bos/alignments/<int:alignment_id>/features/',
    view_func=alignment_feat_view,
    methods=['GET', 'POST']
)
bp_alignments.add_url_rule(
    bcs_api_base + '/bos/alignments/<int:alignment_id>/features/<int:cmt_id>',
    view_func=alignment_feat_view,
    methods=['GET', 'PUT', 'DELETE']
)
