from ..db_models.sysadmin import Identity, SystemFunction, Group, Role, Organization, RoleIdentity
from . import make_simple_rest_crud, parse_request_params

# --------------------------------------------------------------------------------------------------------------------
bp_identities, IdentitiesAPI = make_simple_rest_crud(Identity, "identities")
bp_roles, RolesAPI = make_simple_rest_crud(Role, "roles")
bp_identities_roles, IdentitiesRolesAPI = make_simple_rest_crud(RoleIdentity, "identities_roles")
bp_groups, GroupsAPI = make_simple_rest_crud(Group, "groups")
bp_organizations, OrganizationsAPI = make_simple_rest_crud(Organization, "organizations")
# bp_acl, aclAPI = make_simple_rest_crud(ACL, "acls")

bp_sys_functions, SystemFunctionsAPI = make_simple_rest_crud(SystemFunction, "system_functions")

# --------------------------------------------------------------------------------------------------------------------
# bp_sys_functions = Blueprint('bp_sys_functions', __name__)
#
#
# class SystemFunctionsAPI(RestfulView):
#     primary_key = ('id',)
#
#     @bcs_session(read_only=True)
#     def list(self):
#         db = g.bcs_session.db_session
#         r = ResponseObject()
#         r.content = db.query(SystemFunction).all()
#         return r.get_response()
#
#     @bcs_session()
#     def create(self):
#         db = g.bcs_session.db_session
#         r = ResponseObject()
#         t = request.json
#         s = SystemFunction.Schema().load(t, instance=SystemFunction())
#         db.add(s)
#         return r.get_response()
#
#     @bcs_session()
#     def get(self, _id):
#         db = g.bcs_session.db_session
#         r = ResponseObject()
#         r.content = db.query(SystemFunction).get(_id)
#         return r.get_response()
#
#     @bcs_session()
#     def replace(self, _id):
#         db = g.bcs_session.db_session
#         r = ResponseObject()
#         s = SystemFunction()
#         db.save(s, _id)
#         return r.get_response()
#
#     @bcs_session()
#     def update(self, _id):
#         db = g.bcs_session.db_session
#         r = ResponseObject()
#         s = SystemFunction()
#         db.save(s)
#         return r.get_response()
#
#     @bcs_session()
#     def delete(self, _id):
#         db = g.bcs_session.db_session
#         r = ResponseObject()
#         s = SystemFunction()
#         db.delete(s)
#         return r.get_response()
#
#
# SystemFunctionsAPI.\
#     route_as_view(bp_sys_functions, 'system_functions',
#                   (f"{bcs_api_base}/system_functions/", f"{bcs_api_base}/system_functions/<int:_id>"))

from flask import Blueprint, request
from flask.views import MethodView
from biobarcoding.authentication import n_session
from biobarcoding.rest import app_api_base, ResponseObject
from biobarcoding.services.acls import *

bp_acl = Blueprint('bp_acl', __name__)


class ACLAPI(MethodView):
    """
    ACLs Resource
    """

    @n_session(read_only=True)
    def get(self, id=None, uuid=None):
        print(f'GET {request.path}\nGetting ACL {id}')
        kwargs = parse_request_params()
        issues, content, count, status = read_acls(id=id, uuid=uuid, **kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()

    @n_session()
    def post(self):
        print(f'POST {request.path}\nCreating ACL')
        kwargs = parse_request_params()
        issues, content, status = create_acls(**kwargs.get('values'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    @n_session()
    def put(self, id=None, uuid=None):
        print(f'PUT {request.path}\nCreating ACL {id}')
        kwargs = parse_request_params()
        issues, content, status = update_acls(id=id, uuid=uuid, **kwargs.get('values'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    @n_session()
    def delete(self, id=None, uuid=None):
        print(f'DELETE {request.path}\nDeleting ACL {id}')
        kwargs = parse_request_params()
        issues, content, status = delete_acls(id=id, uuid=uuid, **kwargs)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


acl_view = ACLAPI.as_view('api_acl')
bp_acl.add_url_rule(
    app_api_base + '/acls',
    view_func=acl_view,
    methods=['GET', 'POST', 'DELETE']
)
bp_acl.add_url_rule(
    app_api_base + '/acls/',
    view_func=acl_view,
    methods=['GET', 'POST', 'DELETE']
)
bp_acl.add_url_rule(
    app_api_base + '/acls/<int:id>',
    view_func=acl_view,
    methods=['GET', 'PUT', 'DELETE']
)
bp_acl.add_url_rule(
    app_api_base + '/acls/<string:uuid>',
    view_func=acl_view,
    methods=['GET', 'PUT', 'DELETE']
)


class ObjectTypeAPI(MethodView):
    """
    Object Types Resource
    """

    @n_session(read_only=True)
    def get(self, id=None, uuid=None):
        print(f'GET {request.path}\nGetting ACL {id}')
        kwargs = parse_request_params()
        issues, content, count, status = read_obj_types(id=id, uuid=uuid, **kwargs)
        return ResponseObject(content=content, count=count, issues=issues, status=status).get_response()


obj_type_view = ObjectTypeAPI.as_view('api_obj_type')
bp_acl.add_url_rule(
    app_api_base + '/object_types',
    view_func=obj_type_view,
    methods=['GET']
)
bp_acl.add_url_rule(
    app_api_base + '/object_types/',
    view_func=obj_type_view,
    methods=['GET']
)
bp_acl.add_url_rule(
    app_api_base + '/object_types/<int:id>',
    view_func=obj_type_view,
    methods=['GET']
)
bp_acl.add_url_rule(
    app_api_base + '/object_types/<string:uuid>',
    view_func=obj_type_view,
    methods=['GET']
)
