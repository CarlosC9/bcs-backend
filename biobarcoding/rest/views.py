from flask import Blueprint, request
from flask.views import MethodView

from ..authentication import n_session
from . import app_api_base, ResponseObject
from ..services.views import *

bp_views = Blueprint('bp_views', __name__)


class ViewAPI(MethodView):
    """
    Views Resource
    """

    @n_session(read_only=True)
    def get(self, view=None):
        print(f'GET {request.path}\nGetting {view}')
        kwargs = self._prepare(view)
        issues, content, count, status = read(view, **kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    def _prepare(self, view):
        from biobarcoding.rest import check_request_params
        return check_request_params()


views_view = ViewAPI.as_view('api_views')
bp_views.add_url_rule(
    app_api_base + '/views',
    view_func=views_view,
    methods=['GET']
)
bp_views.add_url_rule(
    app_api_base + '/views/<string:view>',
    view_func=views_view,
    methods=['GET']
)
