import importlib
import os.path

from flask import Blueprint, request, send_file
from flask.views import MethodView

from ..authentication import n_session
from . import app_api_base, ResponseObject, Issue, IType, parse_request_params

bp_bos = Blueprint('bp_bos', __name__)


class BioObjAPI(MethodView):
    """
    BOS Resource
    """

    @n_session(read_only=True)
    def get(self, bos, id=None, format=None):
        print(f'GET {request.path}\nGetting {bos} {id}')
        kwargs = self._prepare(bos)
        count = 0
        if format:
            issues, content, status = self.service.export(id, format=format, **kwargs)
            if content:
                return send_file(content, mimetype=f'text/{format}'), status
        else:
            issues, content, count, status = self.service.read(id, **kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def post(self, bos, format=None):
        print(f'POST {request.path}\nCreating {bos}')
        kwargs = self._prepare(bos)
        if request.files or kwargs.get('value') and 'filesAPI' in kwargs.get('value'):
            issues, content, status = self._import_files(format, kwargs.get('value'))
        else:
            issues, content, status = self.service.create(**kwargs.get('value'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    @n_session()
    def put(self, bos, id, format=None):
        print(f'PUT {request.path}\nUpdating {bos} {id}')
        kwargs = self._prepare(bos)
        issues, content, status = self.service.update(id, **kwargs.get('value'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    @n_session()
    def delete(self, bos, id=None, format=None):
        print(f'DELETE {request.path}\nDeleting {bos} {id}')
        kwargs = self._prepare(bos)
        issues, content, status = self.service.delete(id, **kwargs)
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    def _import_files(self, format, value={}):
        issues, content = [], []
        if 'filesAPI' in value:
            i, c, s = self._import_filesAPI(format, value)
            issues += i
            content += c
        if request.files:
            i, c, s = self._import_request_files(format, value)
            issues += i
            content += c
        return issues, content, 207

    def _import_filesAPI(self, format, value={}):
        issues, content = [], []
        if value.get('filesAPI'):
            # TODO: deal with list of filepaths
            filesAPI = value.get('filesAPI') \
                if isinstance(value.get('filesAPI'), (list, tuple)) \
                else [value.get('filesAPI')]
            from biobarcoding.db_models import DBSession
            from biobarcoding.db_models.files import FileSystemObject
            for file in filesAPI:
                try:
                    if not os.path.isabs(file):
                        raise Exception('Invalid path')
                    file = DBSession.query(FileSystemObject) \
                        .filter(FileSystemObject.full_name == file).first()
                    from werkzeug.utils import secure_filename
                    file_cp = '/tmp/' + secure_filename(file.full_name)
                    with open(file_cp, 'wb') as f:
                        f.write(file.embedded_content)
                    i, c, s = self.service.import_file(file_cp, format, **value)
                except Exception as e:
                    print(e)
                    i, c = [Issue(IType.ERROR, f'Could not import the file {file}.', file)], {}
                issues += i
                content.append(c)
        return issues, content, 207

    def _import_request_files(self, format, value={}):
        issues, content = [], []
        for key, file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                i, c, s = self.service.import_file(file_cpy, format, **value)
            except Exception as e:
                print(e)
                i, c = [Issue(IType.ERROR, f'Could not import the file {file}.', file.filename)], {}
            issues += i
            content.append(c)
        return issues, content, 207

    def _make_file(self, file):
        import os
        from werkzeug.utils import secure_filename
        file_path = os.path.join('/tmp', secure_filename(file.filename))
        file.save(file_path)
        return file_path

    def _prepare(self, bos):
        self.service = importlib.import_module(f'biobarcoding.services.{bos}')
        return parse_request_params()


bos_view = BioObjAPI.as_view('api_bos')
bp_bos.add_url_rule(
    app_api_base + '/bos/<string:bos>',
    view_func=bos_view,
    methods=['GET', 'POST', 'DELETE']
)
bp_bos.add_url_rule(
    app_api_base + '/bos/<string:bos>/',
    view_func=bos_view,
    methods=['GET', 'POST', 'DELETE']
)
bp_bos.add_url_rule(
    app_api_base + '/bos/<string:bos>/<string:id>',
    view_func=bos_view,
    methods=['GET', 'PUT', 'DELETE']
)
bp_bos.add_url_rule(
    app_api_base + '/bos/<string:bos>.<string:format>',
    view_func=bos_view,
    methods=['GET', 'POST']
)
bp_bos.add_url_rule(
    app_api_base + '/bos/<string:bos>/.<string:format>',
    view_func=bos_view,
    methods=['GET', 'POST']
)
bp_bos.add_url_rule(
    app_api_base + '/bos/<string:bos>/<string:id>.<string:format>',
    view_func=bos_view,
    methods=['GET', 'PUT', 'DELETE']
)

# class BioFeatAPI(MethodView):
#     """
#     Sequences Feature Resource
#     """
#     def get(self, seq_id, cmt_id=None):
#         print(f'GET {request.path}\nGetting comments {seq_id} {cmt_id}')
#         return ResponseObject(content={'status':'success','message':'READ: sequence comments, dummy complete'}, status=200).get_response()
#
#
#     def post(self, seq_id):
#         print(f'POST {request.path}\nCreating comments')
#         return ResponseObject(content={'status':'success','message':'CREATE: sequence comments, dummy complete'}, status=200).get_response()
#
#
#     def put(self, seq_id, cmt_id):
#         print(f'PUT {request.path}\nCreating comments {seq_id} {cmt_id}')
#         return ResponseObject(content={'status':'success','message':'UPDATE: sequence comments, dummy complete'}, status=200).get_response()
#
#
#     def delete(self, seq_id, cmt_id=None):
#         print(f'DELETE {request.path}\nDeleting comments {seq_id} {cmt_id}')
#         return ResponseObject(content={'status':'success','message':'DELETE: sequence comments, dummy complete'}, status=200).get_response()
#
#
#     def _check_data(self, data):
#         if data:
#             if data.get('id'):
#                 self.ids = data['id']
#             if data.get('comment'):
#                 self.comment = data['comment']
#             if data.get('range'):
#                 self.range = data['range']
#         print(f'DATA: {data}')
#
#
# seq_feat_view = SeqFeatAPI.as_view('api_seq_feat')
# bp_bos.add_url_rule(
#     bcs_api_base + '/bos/<string:bos>/<int:seq_id>/features/',
#     view_func=seq_feat_view,
#     methods=['GET','POST']
# )
# bp_bos.add_url_rule(
#     bcs_api_base + '/bos/<string:bos>/<int:seq_id>/features/<int:cmt_id>',
#     view_func=seq_feat_view,
#     methods=['GET','PUT','DELETE']
# )
