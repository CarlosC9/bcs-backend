from biobarcoding.db_models.sysadmin import Identity, SystemFunction, Group, Role, Organization, ACL, RoleIdentity
from biobarcoding.rest import make_simple_rest_crud

# --------------------------------------------------------------------------------------------------------------------
bp_identities, IdentitiesAPI = make_simple_rest_crud(Identity, "identities")
bp_roles, RolesAPI = make_simple_rest_crud(Role, "roles")
bp_identities_roles, IdentitiesRolesAPI = make_simple_rest_crud(RoleIdentity, "identities_roles")
bp_groups, GroupsAPI = make_simple_rest_crud(Group, "groups")
bp_organizations, OrganizationsAPI = make_simple_rest_crud(Organization, "organizations")
bp_acl, aclAPI = make_simple_rest_crud(ACL, "acls")

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
