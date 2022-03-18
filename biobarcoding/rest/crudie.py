from flask import request, send_file
from flask.views import MethodView

from . import parse_request_params, ResponseObject, Issue, IType
from ..authentication import n_session
from ..services.main import getCRUDIE


class CrudieAPI(MethodView):
    """
    CRUDIE API Resource
    """
    params = {}
    service = None

    def __init__(self, entity):
        try:
            self.service = getCRUDIE(entity)
        except Exception as e:
            from flask import abort
            abort(400, f'Bad request: invalid URL ({request.path})')  # TODO: abort properly

    def pre_request(self):
        try:
            self.params = parse_request_params()
        except Exception as e:
            from flask import abort
            abort(400, f'Bad request: invalid params')  # TODO: abort properly

    # CRUDIE REQUESTS

    @n_session(read_only=True)
    def get(self, endpoint=None, id=None, format=None, **kwargs):
        print(f'GET {request.path}\nGetting {endpoint}')
        self.pre_request()

        if format:
            issues, content, count, status = self.service.export_file(id=id, format=format, **self.params)
            if content:
                return send_file(content, mimetype=f'text/{format}'), status
        else:
            issues, content, count, status = self.service.read(id=id, format=format, **self.params)

        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def post(self, endpoint=None, **kwargs):
        print(f'POST {request.path}\nCreating {endpoint}')
        self.pre_request()

        if request.files or self.params.get('values', {}).get('filesAPI'):
            issues, content, count, status = self._import_files(**self.params.get('values'))
        else:
            issues, content, count, status = self.service.create(**self.params.get('values'))

        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def put(self, endpoint=None, id=None, **kwargs):
        print(f'PUT {request.path}\nUpdating {endpoint} {id}')
        self.pre_request()

        issues, content, count, status = self.service.update(id=id, **self.params)

        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def delete(self, endpoint=None, id=None, **kwargs):
        print(f'DELETE {request.path}\nDeleting {endpoint}')
        self.pre_request()

        issues, content, count, status = self.service.delete(id=id, **self.params)

        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    # ADDITIONAL TOOLS

    def _import_files(self, **values):
        # TODO: allow combinations ?
        if values.get('filesAPI'):
            return self._import_filesAPI(**values)
        if request.files:
            return self._import_request_files(**values)
        return [], None, 0, 200

    def _import_filesAPI(self, **values):
        issues, content = [], []
        if values.get('filesAPI'):
            file = values.get('filesAPI')
            from biobarcoding.db_models import DBSession
            from biobarcoding.db_models.files import FileSystemObject
            try:
                import os
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
                issues, content, count, status = self.service.import_file(file_cp, **values)
            except Exception as e:
                print(e)
                issues += [Issue(IType.ERROR, f'Could not import the file {file}.', file)]
        return issues, content, count, 207

    def _import_request_files(self, **values):
        issues, content, count = [], [], 0
        for key, file in request.files.items(multi=True):
            try:
                file_cpy = self._make_file(file)
                i, c, s = self.service.import_file(file_cpy, **values)
                count += 1
            except Exception as e:
                print(e)
                i, c = [Issue(IType.ERROR, f'Could not import the file {file}.', file.filename)], {}
            issues += i
            content.append(c)
        return issues, content, count, 207

    def _make_file(self, file):
        import os
        from werkzeug.utils import secure_filename
        file_path = os.path.join('/tmp', secure_filename(file.filename))
        file.save(file_path)
        return file_path

