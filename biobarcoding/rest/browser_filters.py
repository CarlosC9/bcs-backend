from flask import Blueprint, request
from flask.views import MethodView

from ..authentication import n_session
from . import bcs_api_base, ResponseObject, check_request_params
from ..services import browser_filters

bp_bfilters = Blueprint('bp_bfilters', __name__)


class BrowserFilterAPI(MethodView):
    """
    Browser filters Resource
    """
    kwargs = {}

    @n_session(read_only=True)
    def get(self, datatype, id=None):
        print(f'GET {request.path}\nGetting filters {id}')
        self.kwargs = check_request_params()
        issues, content, count, status = browser_filters.read(datatype, id, **self.kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def post(self, datatype):
        print(f'POST {request.path}\nCreating filters')
        self.kwargs = check_request_params()
        issues, content, status = browser_filters.create(datatype, **self.kwargs.get('value'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    @n_session()
    def put(self, datatype, id):
        print(f'PUT {request.path}\nUpdating filter {id}')
        self.kwargs = check_request_params()
        issues, content, status = browser_filters.update(datatype, id, **self.kwargs.get('value'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    @n_session()
    def delete(self, datatype, id=None):
        print(f'DELETE {request.path}\nDeleting filters {id}')
        self.kwargs = check_request_params()
        issues, content, status = browser_filters.delete(datatype, id, **self.kwargs)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


bfilter_view = BrowserFilterAPI.as_view('api_bfilter')
bp_bfilters.add_url_rule(
    bcs_api_base + '/browser/<string:datatype>/filters/',
    view_func=bfilter_view,
    methods=['GET', 'POST', 'DELETE']
)
bp_bfilters.add_url_rule(
    bcs_api_base + '/browser/<string:datatype>/filters/<string:id>',
    view_func=bfilter_view,
    methods=['GET', 'PUT', 'DELETE']
)


class BrowserFilterFormAPI(MethodView):
    """
    Browser filter forms Resource
    """

    def get(self, datatype):
        print(f'GET {request.path}\nGetting filter forms {datatype}')
        issues, content, status = browser_filters.read_form(datatype)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


bfilter_forms_view = BrowserFilterFormAPI.as_view('api_bfilter_forms')
bp_bfilters.add_url_rule(
    bcs_api_base + '/browser/<string:datatype>/filters/schema',
    view_func=bfilter_forms_view,
    methods=['GET']
)
