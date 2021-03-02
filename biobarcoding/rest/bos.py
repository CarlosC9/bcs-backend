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
        bos_service = self._prepare(bos)
        if format:
            issues, content, status = bos_service.export(id, format=format, **self.kwargs)
            return send_file(content, mimetype=f'text/{format}'), status
        issues, content, status = bos_service.read(id, **self.kwargs)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def post(self, bos, format=None):
        print(f'POST {request.path}\nCreating {bos}')
        if request.files:
            issues, content, status = self._import_files(format)
        else:
            bos_service = self._prepare(bos)
            issues, content, status = bos_service.create(**self.kwargs)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def put(self, bos, id, format=None):
        print(f'PUT {request.path}\nCreating {bos} {id}')
        bos_service = self._prepare(bos)
        issues, content, status = bos_service.update(id, **self.kwargs)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    @bcs_session()
    def delete(self, bos, id=None, format=None):
        print(f'DELETE {request.path}\nDeleting {bos} {id}')
        bos_service = self._prepare(bos)
        issues, content, status = bos_service.delete(id, **self.kwargs)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


    def _import_files(self, bos, format):
        bos_service = self._prepare(bos)
        issues, content = [], []
        for key,file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                i, c, s = bos_service.import_file(file_cpy, format, **self.kwargs)
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
        return importlib.import_module(f'biobarcoding.services.{bos}')


    def _check_data(self, data=None):
        if not data:
            self.kwargs = { 'filter' : [],
              'order' : {},
              'pagination' : {},
              'value' : {}}
            if request.json:
                self._check_data(request.json)
            if request.values:
                self._check_data(request.values)
        else:
            print(f'DATA: {data}')
            # dict = data.to_dict(flat=False)
            # args = { key: dict[key][0] if len(dict[key]) <= 1 else dict[key] for key in dict }
            # self.kwargs.update(args)
            if 'filter' in data and data['filter']:
                from urllib.parse import unquote
                self.kwargs['filter'] += json.loads(unquote(data.get('filter')))
            if 'order' in data and data['order']:
                self.kwargs['order'].update(
                    json.loads(
                        data.get('order')))
            if 'pagination' in data and data['pagination']:
                self.kwargs['pagination'].update(
                    json.loads(
                        data.get('pagination')))
            if 'value' in data and data['value']:
                self.kwargs['value'].update(
                    json.loads(
                        data.get('value')))
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
#             if 'id' in data and data['id']:
#                 self.ids = data['id']
#             if 'comment' in data and data['comment']:
#                 self.comment = data['comment']
#             if 'range' in data and data['range']:
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
