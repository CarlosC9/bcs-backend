from ..db_models.sysadmin import Identity, SystemFunction, Group, Role, Organization, RoleIdentity, \
    IdentityAuthenticator
from . import make_simple_rest_crud, parse_request_params, tm_default_users

from flask import Blueprint, request
from flask.views import MethodView
from ..authentication import n_session
from ..rest import app_api_base, ResponseObject
from ..services.acls import *


# --------------------------------------------------------------------------------------------------------------------

def custom_identities_filter(filter, session=None):
    """
    Obtain non-system identities

    :param filter: A dictionary (str, dict)
    :return:
    """
    return [Identity.uuid.notin_([i for i in tm_default_users.keys()])]


bp_identities, IdentitiesAPI = make_simple_rest_crud(Identity, "identities", aux_filter=custom_identities_filter, default_filter={'-': {}})
bp_identities_authenticators, IdentitiesAuthenticatorAPI = make_simple_rest_crud(IdentityAuthenticator, "identities_authenticators")
bp_roles, RolesAPI = make_simple_rest_crud(Role, "roles")
bp_identities_roles, IdentitiesRolesAPI = make_simple_rest_crud(RoleIdentity, "identities_roles")
bp_groups, GroupsAPI = make_simple_rest_crud(Group, "groups")
bp_organizations, OrganizationsAPI = make_simple_rest_crud(Organization, "organizations")
# bp_acl, aclAPI = make_simple_rest_crud(ACL, "acls")

bp_sys_functions, SystemFunctionsAPI = make_simple_rest_crud(SystemFunction, "system_functions")

# --------------------------------------------------------------------------------------------------------------------

bp_acl = Blueprint('bp_acl', __name__)


class ACLAPI(MethodView):
    """
    ACLs Resource
    """

    @n_session(read_only=True)
    def get(self, id=None, uuid=None):
        print(f'GET {request.path}\nGetting ACL {id}')
        kwargs = parse_request_params()
        issues, content, count, status = read_acls(id_=id, uuid=uuid, **kwargs)
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
        issues, content, status = update_acls(id_=id, uuid=uuid, **kwargs.get('values'))
        return ResponseObject(content=content, issues=issues, status=status).get_response()

    @n_session()
    def delete(self, id=None, uuid=None):
        print(f'DELETE {request.path}\nDeleting ACL {id}')
        kwargs = parse_request_params()
        issues, content, status = delete_acls(id_=id, uuid=uuid, **kwargs)
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
        print(f'GET {request.path}\nGetting Object Types {id}')
        kwargs = parse_request_params()
        issues, content, count, status = read_obj_types(id_=id, uuid=uuid, **kwargs)
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
