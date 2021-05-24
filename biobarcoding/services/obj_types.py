from biobarcoding.db_models import DBSession, ObjectType
from biobarcoding.rest import Issue, IType
from biobarcoding.services import get_query


def create(**kwargs):
    issues, status = [Issue(IType.INFO, f'CREATE obj_types: dummy successfully completed.')], 200
    return issues, None, 0, status


count = 0
def read(id=None, type=None, **kwargs):
    content = None
    try:
        content = get_query(DBSession, ObjectType, id=id, name=type, **kwargs)
        if id:
            content = content.first()
        elif type:
            content = __get_permission(content).first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ obj_types: The obj_types were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ obj_types: The obj_types could not be read.')], 500
    return issues, content, count, status


def update(id, type=None, **kwargs):
    issues, status = [Issue(IType.INFO, f'UPDATE obj_types({type}): dummy successfully completed.')], 200
    return issues, None, status


def delete(id=None, type=None, **kwargs):
    issues, status = [Issue(IType.INFO, f'DELETE obj_types({type}): dummy successfully completed.')], 200
    return issues, None, status


def __get_permission(query):
    # ObjectType.name > ACL.object_type > ACLDetail.acl_id > ACLDetail.permission_id > PermissionType
    from biobarcoding.db_models.sysadmin import ACL, ACLDetail, PermissionType
    from sqlalchemy.orm import load_only
    return DBSession.query(ACLDetail.permission).filter(
        ACLDetail.acl.has(ACL.object_type.in_(query.options(load_only("id")))))
    # return DBSession.query(PermissionType).filter(PermissionType.id.in_(ids))
