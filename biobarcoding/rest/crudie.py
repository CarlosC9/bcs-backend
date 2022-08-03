import os.path

from flask import request, send_file, current_app
from flask.views import MethodView
from werkzeug.utils import secure_filename

from . import parse_request_params, ResponseObject, Issue, IType
from ..authentication import n_session
from ..services import log_exception
from ..services.main import getCRUDIE


class CrudieAPI(MethodView):
    """
    CRUDIE API Resource
    """
    params = {}
    service = None

    def __init__(self, entity):
        try:
            self.entity = entity
            self.service = getCRUDIE(self.entity)
        except Exception as e:
            print(f'Bad request: invalid URL ({request.path})')
            log_exception(e)
            from flask import abort
            abort(400, f'Bad request: invalid URL ({request.path})')  # TODO: abort properly

    def pre_request(self, **kwargs):
        try:
            self.params = parse_request_params()
            self.params.update(kwargs)
        except Exception as e:
            print(f'Bad request: invalid params')
            log_exception(e)
            from flask import abort
            abort(400, f'Bad request: invalid params')  # TODO: abort properly

    # CRUDIE REQUESTS

    @n_session(read_only=True)
    def get(self, id=None, format=None, **kwargs):
        print(f'GET {request.path}\nGetting {self.entity}')
        self.pre_request(id=id, format=format, **kwargs)

        if format:
            issues, content, count, status = self.service.export_data(**self.params)
            if content:
                return send_file(content, mimetype=f'text/{format}'), status
        else:
            issues, content, count, status = self.service.read(**self.params)

        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def post(self, **kwargs):
        print(f'POST {request.path}\nCreating {self.entity}')
        self.pre_request(**kwargs)

        if self.params.get('values', {}).get('filesAPI'):
            issues, content, count, status = self._import_filesAPI(**self.params.get('values'))
        elif request.files:
            issues, content, count, status = self._import_request_files(**self.params.get('values'))
        else:
            issues, content, count, status = self.service.create(**self.params.get('values'))

        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def put(self, id=None, **kwargs):
        print(f'PUT {request.path}\nUpdating {self.entity} {id}')
        self.pre_request(id=id, **kwargs)

        # TODO update from files
        issues, content, count, status = self.service.update(**self.params)

        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def delete(self, id=None, **kwargs):
        print(f'DELETE {request.path}\nDeleting {self.entity}')
        self.pre_request(id=id, **kwargs)

        issues, content, count, status = self.service.delete(**self.params)

        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    # ADDITIONAL TOOLS

    def _prepare_filesAPI(self, file):
        from ..db_models import DBSession
        from ..db_models.files import FileSystemObject
        if not os.path.isabs(file):
            raise Exception('Invalid filepath')
        # TODO: might full_name not be unique
        file = DBSession.query(FileSystemObject) \
            .filter(FileSystemObject.full_name == file).first()
        file_cp = os.path.join(current_app.config['UPLOAD_FOLDER'], secure_filename(file.full_name))
        with open(file_cp, 'wb') as f:
            f.write(file.embedded_content)
        return file_cp

    def _import_filesAPI(self, **values):
        issues, content, count = [], [], 0
        files = values.pop('filesAPI', None)
        if not isinstance(files, (tuple, list)):
            files = [files]
        for file in files:
            try:
                file_cp = self._prepare_filesAPI(file)
                i, c, cc, s = self.service.import_data(file_cp, filesAPI=file, **values)
            except Exception as e:
                print(e)
                i, c, cc = [Issue(IType.ERROR, f'Could not import the file {file}.', file)], {}, 0
            issues += i
            content.append(c)
            count += cc
        return issues, content, count, 207

    def _import_request_files(self, **values):
        issues, content, count = [], [], 0
        for key, file in request.files.items(multi=True):
            try:
                file_cp = os.path.join(current_app.config['UPLOAD_FOLDER'], secure_filename(file.filename))
                file.save(file_cp)
                i, c, cc, s = self.service.import_data(file_cp, **values)
            except Exception as e:
                print(e)
                i, c, cc = [Issue(IType.ERROR, f'Could not import the file {file}.', file.filename)], {}, 0
            issues += i
            content.append(c)
            count += cc
        return issues, content, count, 207

