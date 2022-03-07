from sqlalchemy import and_

from ..db_models.sysadmin import Identity, SystemFunction, Group, Role, Organization, RoleIdentity, \
    IdentityAuthenticator, IdentityStoreEntry
from . import make_simple_rest_crud, parse_request_params, tm_default_users

from flask import Blueprint, request, g
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


def get_current_identity(session, identity_id):
    if identity_id is not None:
        if int(identity_id) == 0:
            identity_id = g.n_session.identity.id
            return session.query(Identity).filter(Identity.id == identity_id).first()


bp_identities, IdentitiesAPI = make_simple_rest_crud(Identity, "identities", alt_getters=dict(get=get_current_identity),
                                                     aux_filter=custom_identities_filter, default_filter={'-': {}})
bp_identities_authenticators, IdentitiesAuthenticatorAPI = make_simple_rest_crud(IdentityAuthenticator, "identities_authenticators")
bp_roles, RolesAPI = make_simple_rest_crud(Role, "roles")


def get_identity_roles(session, identity_id):
    """
    Load roles for a given (or not) identity

    :param session:
    :param identity_id:
    :return:
    """
    if identity_id is not None:
        if int(identity_id) == 0:
            identity_id = g.n_session.identity.id
        return session.query(Role).join(RoleIdentity).filter(RoleIdentity.identity_id == identity_id).all()


bp_identities_roles, IdentitiesRolesAPI = make_simple_rest_crud(RoleIdentity, "identities_roles",
                                                                alt_getters=dict(get=get_identity_roles))
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


class IdentityStoreAPI(MethodView):
    """
    Endpoint to enable storing per identity key-values, where values are stored as JSON

    LOGIN (entries assume the current user):
    export API_BASE_URL=http://localhost:5000/api
    curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
    PUT:
    curl --cookie app-cookies.txt -i -XPUT "http://localhost:5000/api/identity_store/status" -H "Content-Type: application/json" -d '{"status": "scheduled"}'
    GET (all keys):
    curl --cookie app-cookies.txt "http://localhost:5000/api/identity_store/"
    GET (key "status"):
    curl --cookie app-cookies.txt "http://localhost:5000/api/identity_store/status"

    """
    @n_session(read_only=True, authr=None)
    def get(self, key=None):  # List or Read
        db = g.n_session.db_session
        ident = g.n_session.identity
        valid_identity = ident is not None  # and ident != anonymous_identity
        r = ResponseObject()
        if not valid_identity:
            r.status = 401
            r.issues.append('Not identified, login first')
        else:
            if key is None:
                # List of all
                kwargs = parse_request_params()
                from ..services import get_query
                # Filter always by identity (override possible "filter", not "order" or others)
                kwargs.update(dict(filter=[dict(identity_id={'op': 'eq', 'unary': ident.id})]))
                query, count = get_query(db, IdentityStoreEntry, **kwargs)
                r.count = count
                r.content = [e.key for e in query.all()]
                # r.content = {e.key: e.value for e in query.all()}
            else:
                # Detail
                _ = db.query(IdentityStoreEntry).\
                    filter(and_(IdentityStoreEntry.identity == ident, IdentityStoreEntry.key == key)).first()
                if _ is None:
                    r.status = 404
                    r.issues.append(Issue(IType.ERROR, f'Key {key} not found'))
                else:
                    r.content = _.value
                    r.count = 1

        return r.get_response()

    @n_session()
    def put(self, key):  # Update (total or partial)
        db = g.n_session.db_session
        ident = g.n_session.identity
        valid_identity = ident is not None  # and ident != anonymous_identity
        r = ResponseObject()
        if not valid_identity:
            r.status = 401
            r.issues.append('Not identified, login first')
        else:
            s = db.query(IdentityStoreEntry).\
                filter(and_(IdentityStoreEntry.identity == ident, IdentityStoreEntry.key == key)).first()
            if s is None:
                s = IdentityStoreEntry()
                db.add(s)
                s.identity = ident
                s.key = key
            s.value = request.json
            r.content = s
        return r.get_response()

    @n_session()
    def delete(self, key):  # Delete
        db = g.n_session.db_session
        ident = g.n_session.identity
        valid_identity = ident is not None  # and ident != anonymous_identity
        r = ResponseObject()
        if not valid_identity:
            r.status = 401
            r.issues.append('Not identified, login first')
        else:
            s = db.query(IdentityStoreEntry).\
                filter(and_(IdentityStoreEntry.identity == ident, IdentityStoreEntry.key == key)).first()
            db.delete(s)
        return r.get_response()


bp_identity_store = Blueprint(f'bp_identity_store', __name__)
view_func = IdentityStoreAPI.as_view("identity_store")
bp_identity_store.add_url_rule(f"{app_api_base}/identity_store/<string:key>", view_func=view_func, methods=['GET', 'PUT', 'DELETE'])
bp_identity_store.add_url_rule(f"{app_api_base}/identity_store/", defaults=dict(key=None), view_func=view_func, methods=['GET'])