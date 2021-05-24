from flask import Blueprint, request
from flask.views import MethodView
from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject

from biobarcoding.services.obj_types import *

bp_obj_types = Blueprint('bp_obj_types', __name__)


class ObjTypeAPI(MethodView):
    """
    Object Types Resource
    """
    @bcs_session(read_only=True)
    def get(self, type=None):
        print(f'GET {request.path}\nGetting {type}')
        kwargs = self._prepare(type)
        issues, content, count, status = read(type=type, **kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()


    def _prepare(self, type):
        from biobarcoding.rest import check_request_params
        return check_request_params()


obj_types_view = ObjTypeAPI.as_view('api_obj_types')
bp_obj_types.add_url_rule(
    bcs_api_base + '/object_types',
    view_func=obj_types_view,
    methods=['GET']
)
bp_obj_types.add_url_rule(
    bcs_api_base + '/object_types/<string:type>',
    view_func=obj_types_view,
    methods=['GET']
)