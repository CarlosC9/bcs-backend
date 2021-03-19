import importlib, json
from flask import Blueprint, request, send_file
from flask.views import MethodView
from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject, Issue, IType


bp_bos = Blueprint('bp_bos', __name__)


class BioObjAPI(MethodView):
    """
    BOS Resource
    """
    kwargs = {}


    @bcs_session(read_only=True)
    def get(self, bos, id=None, format=None):
        print(f'GET {request.path}\nGetting {bos} {id}')
        self._prepare(bos)
        if format:
            issues, content, status = self.service.export(id, format=format, **self.kwargs)
            return send_file(content, mimetype=f'text/{format}'), status
        issues, content, count, status = self.service.read(id, **self.kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()


    @bcs_session()
    def post(self, bos, format=None):
        print(f'POST {request.path}\nCreating {bos}')
        self._prepare(bos)
        if request.files:
            issues, content, status = self._import_files(format)
        else:
            issues, content, status = self.service.create(**self.kwargs.get('value'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def put(self, bos, id, format=None):
        print(f'PUT {request.path}\nCreating {bos} {id}')
        self._prepare(bos)
        issues, content, status = self.service.update(id, **self.kwargs.get('value'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def delete(self, bos, id=None, format=None):
        print(f'DELETE {request.path}\nDeleting {bos} {id}')
        self._prepare(bos)
        issues, content, status = self.service.delete(id, **self.kwargs)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    def _import_files(self, format):
        issues, content = [], []
        for key,file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                i, c, s = self.service.import_file(file_cpy, format, **self.kwargs.get('value'))
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


    def _prepare(self, bos):
        self._check_data()
        self.service = importlib.import_module(f'biobarcoding.services.{bos}')


    def _check_data(self, data=None):
        if not data:
            self.kwargs = { 'filter' : [], 'order' : {}, 'pagination' : {},
                            'value' : {}, 'searchValue' : '' }
            if request.json:
                self._check_data(request.json)
            if request.values:
                self._check_data(request.values)
        else:
            print(f'DATA: {data}')
            input = data.copy()
            from urllib.parse import unquote
            for key in self.kwargs:
                if input.get(key):
                    i = input.get(key)
                    try:
                        i = unquote(i)
                    except Exception as e:
                        pass
                    try:
                        i = json.loads(i)
                    except Exception as e:
                        pass
                    self.kwargs[key] = i
                    input.pop(key)
            self.kwargs['value'].update(input)
            print(f'KWARGS: {self.kwargs}')


bos_view = BioObjAPI.as_view('api_bos')
bp_bos.add_url_rule(
    bcs_api_base + '/bos/<string:bos>/',
    view_func=bos_view,
    methods=['GET','POST','DELETE']
)
bp_bos.add_url_rule(
    bcs_api_base + '/bos/<string:bos>/<string:id>',
    view_func=bos_view,
    methods=['GET','PUT','DELETE']
)
bp_bos.add_url_rule(
    bcs_api_base + '/bos/<string:bos>.<string:format>',
    view_func=bos_view,
    methods=['GET','POST']
)
bp_bos.add_url_rule(
    bcs_api_base + '/bos/<string:bos>/<string:id>.<string:format>',
    view_func=bos_view,
    methods=['GET','PUT','DELETE']
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
