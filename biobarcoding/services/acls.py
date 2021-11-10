from ..db_models import DBSession, ObjectType
from ..db_models.bioinformatics import BioinformaticObject
from ..db_models.sysadmin import ACL, ACLDetail, PermissionType
from ..rest import Issue, IType
from . import orm2json, get_simple_query, get_query


def create_acls(**kwargs):
    content = None
    try:
        if kwargs.get('object_uuid'):
            if not kwargs.get('object_type'):
                # TODO: what for the non-bos ?
                kwargs['object_type'] = DBSession.query(BioinformaticObject.bo_type_id) \
                    .filter(BioinformaticObject.uuid == kwargs.get('object_uuid')).first()
        elif kwargs.get('object_type'):
            if not kwargs.get('chado_id'):
                raise Exception('Missing the chado_id')
            kwargs['object_uuid'] = DBSession.query(BioinformaticObject) \
                .filter(BioinformaticObject.chado_id == kwargs.pop('chado_id'),
                        BioinformaticObject.bo_type_id == kwargs.get('object_type')).one().uuid
        else:
            raise Exception('Missing the object_uuid or the chado_id with object_type')
        details = kwargs.pop('details') if 'details' in kwargs else None
        acl = ACL(**kwargs)
        DBSession.add(acl)
        DBSession.flush()
        if isinstance(details, (list, tuple)):
            for d in details:
                DBSession.add(ACLDetail(acl_id=acl.id, **d))
        issues, status = [Issue(IType.INFO, f'CREATE acls: The acl was created successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE acls: The acl could not be created.')], 409
    return issues, content, status


def read_acls(id=None, uuid=None, **kwargs):
    content, count = None, 0
    try:
        qparams = kwargs['values'] if kwargs.get('values') else kwargs
        # IDs: id, uuid, object_uuid, chado_id + object_type
        if qparams.get('chado_id') and not id and not uuid and not qparams.get('object_uuid'):
            if not qparams.get('object_type'):
                raise Exception('Missing the object_type')
            qparams['object_uuid'] = DBSession.query(BioinformaticObject) \
                .filter(BioinformaticObject.chado_id == qparams.pop('chado_id'),
                        BioinformaticObject.bo_type_id == qparams.get('object_type')).one().uuid
        content, count = get_query(DBSession, ACL, id=id, uuid=uuid, **kwargs)

        if id or uuid or qparams.get('object_uuid'):
            content = content.one()
            details = []
            for detail in content.details:
                d = orm2json(detail)
                d['authorizable'] = detail.authorizable
                details.append(d)
            content = orm2json(content)
            content['details'] = details
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ acls: The acls were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ acls: The acls could not be read.')], 404
    return issues, content, count, status


def update_acls(id=None, uuid=None, **kwargs):
    content = None
    try:
        details = kwargs.pop('details')
        content = get_simple_query(DBSession, ACL, id=id, uuid=uuid)
        if kwargs:
            content.update(kwargs)
        content = content.one()
        if isinstance(details, (list, tuple)):
            # Removing missing details
            DBSession.query(ACLDetail).filter(ACLDetail.acl_id == content.id) \
                .filter(ACLDetail.id.notin_([d.get('id') for d in details if d.get('id')])) \
                .delete(synchronize_session='fetch')
            for d in details:
                if 'id' in d:
                    # Updating existing details
                    DBSession.query(ACLDetail).filter(ACLDetail.id == d.get('id')).update(d)
                else:
                    # Adding new details
                    DBSession.add(ACLDetail(acl_id=content.id, **d))
        issues, status = [Issue(IType.INFO, f'UPDATE acls({content.id}): The acls were successfully updated.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'UPDATE acls: The acls could not be updated.')], 409
    return issues, content, status


def delete_acls(id=None, uuid=None, **kwargs):
    content = 0
    try:
        content, count = get_query(DBSession, ACL, id=id, uuid=uuid, **kwargs)
        content = content.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE acls({id}): The acls were successfully deleted.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE acls: The acls could not be deleted.')], 409
    return issues, content, status


"""
ObjectTypes service
"""


def read_obj_types(id=None, uuid=None, **kwargs):
    content, count = None, 0
    try:
        content, count = get_query(DBSession, ObjectType, id=id, uuid=uuid, **kwargs)
        if id or uuid or content.count() == 1:
            content = orm2json(content.first())
            from biobarcoding.db_models.sysadmin import ObjectTypePermissionType
            perms = DBSession.query(PermissionType).join(ObjectTypePermissionType) \
                .filter(ObjectTypePermissionType.object_type_id == content.get('id'))
            content['permission_types'] = perms.all()
            content = [content] if not id and not uuid else content
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ object_types: The object_types were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ object_types: The object_types could not be read.')], 404
    return issues, content, count, status


"""
PermissionTypes service
"""


def read_perm_types(id=None, uuid=None, type=None, **kwargs):
    content, count = None, 0
    try:
        content, count = get_simple_query(DBSession, PermissionType, id=id, uuid=uuid, name=type)
        if id or type:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ permission_types: The permission_types were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ permission_types: The permission_types could not be read.')], 404
    return issues, content, count, status
