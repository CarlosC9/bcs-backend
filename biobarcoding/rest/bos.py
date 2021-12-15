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
        if request.files or kwargs.get('values') and 'filesAPI' in kwargs.get('values'):
            issues, content, status = self._import_files(format, kwargs.get('values'))
        else:
            issues, content, status = self.service.create(**kwargs.get('values'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    @n_session()
    def put(self, bos, id, format=None):
        print(f'PUT {request.path}\nUpdating {bos} {id}')
        kwargs = self._prepare(bos)
        issues, content, status = self.service.update(id, **kwargs.get('values'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    @n_session()
    def delete(self, bos, id=None, format=None):
        print(f'DELETE {request.path}\nDeleting {bos} {id}')
        kwargs = self._prepare(bos)
        issues, content, status = self.service.delete(id, **kwargs)
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    def _import_files(self, format, values={}):
        # TODO: allow combinations ?
        if 'filesAPI' in values:
            return self._import_filesAPI(format, values)
        if request.files:
            return self._import_request_files(format, values)
        return [], None, 200

    def _import_filesAPI(self, format, values={}):
        issues, content = [], []
        if values.get('filesAPI'):
            file = values.get('filesAPI')
            from biobarcoding.db_models import DBSession
            from biobarcoding.db_models.files import FileSystemObject
            try:
                if not os.path.isabs(file):
                    raise Exception('Invalid path')
                if not values.get('sourceuri'):
                    values['sourceuri'] = file
                file = DBSession.query(FileSystemObject) \
                    .filter(FileSystemObject.full_name == file).first()
                from werkzeug.utils import secure_filename
                file_cp = '/tmp/' + secure_filename(file.full_name)
                with open(file_cp, 'wb') as f:
                    f.write(file.embedded_content)
                issues, content, status = self.service.import_file(file_cp, format, **values)
            except Exception as e:
                print(e)
                issues += [Issue(IType.ERROR, f'Could not import the file {file}.', file)]
        return issues, content, 207

    def _import_request_files(self, format, values={}):
        issues, content = [], []
        for key, file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                i, c, s = self.service.import_file(file_cpy, format, **values)
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
        try:
            self.service = importlib.import_module(f'biobarcoding.services.{bos}')
        except Exception as e:
            # TODO: abort properly
            from flask import abort
            abort(400, f'Bad request: unknown {bos} bioinformatic object.')
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
