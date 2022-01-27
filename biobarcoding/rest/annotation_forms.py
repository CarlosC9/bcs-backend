from flask import Blueprint, request
from flask.views import MethodView

from ..authentication import n_session
from . import app_api_base, ResponseObject, parse_request_params
from ..services.main import getCRUDIE

bp_annotations = Blueprint('bp_annotations', __name__)


class AnnotationFormItemAPI(MethodView):
    """
    Annotation form items Resource
    """
    kwargs = {}

    @n_session(read_only=True)
    def get(self, item_type, id=None, cvterm_id=None, db=None, dbxref=None):
        print(f'GET {request.path}\nGetting annotation {item_type}')
        self.kwargs = parse_request_params()
        issues, content, count, status = getCRUDIE(item_type).read(id=id, cvterm_id=cvterm_id, db=db, dbxref=dbxref, **self.kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def post(self, item_type):
        print(f'POST {request.path}\nCreating annotation {item_type}')
        self.kwargs = parse_request_params()
        issues, content, count, status = getCRUDIE(item_type).create(**self.kwargs.get('values'))
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def put(self, item_type, id=None, cvterm_id=None, db=None, dbxref=None):
        print(f'PUT {request.path}\nUpdating annotation {item_type}')
        self.kwargs = parse_request_params()
        issues, content, count, status = getCRUDIE(item_type).update(id=id, cvterm_id=cvterm_id, db=db, dbxref=dbxref, **self.kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def delete(self, item_type, id=None, cvterm_id=None, db=None, dbxref=None):
        print(f'DELETE {request.path}\nDeleting annotation {item_type}')
        self.kwargs = parse_request_params()
        issues, content, count, status = getCRUDIE(item_type).delete(id=id, cvterm_id=cvterm_id, db=db, dbxref=dbxref, **self.kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()


annotation_form_item_view = AnnotationFormItemAPI.as_view('api_annotation_form_item')
bp_annotations.add_url_rule(
    app_api_base + '/annotation_form_<string:item_type>/',
    view_func=annotation_form_item_view,
    methods=['GET', 'POST']
)
bp_annotations.add_url_rule(
    app_api_base + '/annotation_form_<string:item_type>/<int:id>',
    view_func=annotation_form_item_view,
    methods=['GET', 'PUT', 'DELETE']
)
bp_annotations.add_url_rule(
    app_api_base + '/annotation_form_<string:item_type>/cvterm:<int:cvterm_id>',
    view_func=annotation_form_item_view,
    methods=['GET', 'PUT', 'DELETE']
)
bp_annotations.add_url_rule(
    app_api_base + '/annotation_form_<string:item_type>/<string:db>:<string:dbxref>',
    view_func=annotation_form_item_view,
    methods=['GET', 'PUT', 'DELETE']
)


class AnnotationItemAPI(MethodView):
    """
    Annotation items Resource
    """
    kwargs = {}

    @n_session(read_only=True)
    def get(self, object_uuid=None, id=None):
        print(f'GET {request.path}\nGetting annotations {id}')
        self.kwargs = parse_request_params()
        issues, content, count, status = getCRUDIE('annotations').read(object_uuid=object_uuid, id=id, **self.kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def post(self, object_uuid=None):
        print(f'POST {request.path}\nCreating annotation')
        self.kwargs = parse_request_params()
        issues, content, count, status = getCRUDIE('annotations').create(object_uuid=object_uuid, **self.kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def put(self, object_uuid=None, id=None):
        print(f'PUT {request.path}\nUpdating annotation {id}')
        self.kwargs = parse_request_params()
        issues, content, count, status = getCRUDIE('annotations').update(object_uuid=object_uuid, id=id, **self.kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def delete(self, object_uuid=None, id=None):
        print(f'DELETE {request.path}\nDeleting annotation {id}')
        self.kwargs = parse_request_params()
        issues, content, count, status = getCRUDIE('annotations').delete(object_uuid=object_uuid, id=id, **self.kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()


annotation_item_view = AnnotationItemAPI.as_view('api_annotation_item')
bp_annotations.add_url_rule(
    app_api_base + '/annotations/',
    view_func=annotation_item_view,
    methods=['GET','POST']
)
bp_annotations.add_url_rule(
    app_api_base + '/annotations/<int:id>',
    view_func=annotation_item_view,
    methods=['GET','PUT','DELETE']
)
bp_annotations.add_url_rule(
    app_api_base + '/annotations/<string:object_uuid>',
    view_func=annotation_item_view,
    methods=['GET','POST','PUT','DELETE']
)
