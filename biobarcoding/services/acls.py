from biobarcoding.db_models import DBSession, ObjectType
from biobarcoding.db_models.bioinformatics import BioinformaticObject
from biobarcoding.db_models.sysadmin import ACL, ACLDetail, PermissionType
from biobarcoding.rest import Issue, IType
from biobarcoding.rest import get_query


def create_acls(**kwargs):
    content = None
    try:
        if kwargs.get('object_id'):
            if not kwargs.get('object_type'):
                # TODO: what for the non-bos ?
                kwargs['object_type'] = DBSession.query(BioinformaticObject.bo_type_id)\
                    .filter(BioinformaticObject.uuid==kwargs.get('object_id')).first()
        elif kwargs.get('object_type'):
            if not kwargs.get('chado_id'):
                raise
            kwargs['object_id'] = DBSession.query(BioinformaticObject)\
                .filter(BioinformaticObject.chado_id==kwargs.pop('chado_id'),
                        BioinformaticObject.bo_type_id==kwargs.get('object_type')).one().uuid
        else:
            raise
        acl = ACL(**kwargs)
        DBSession.add(acl)
        DBSession.flush()
        details = kwargs.get('details')
        if isinstance(details, (list, tuple)):
            for d in details:
                DBSession.add(ACLDetail(acl_id=acl.id, **d))
        elif details:
            DBSession.add(ACLDetail(acl_id=acl.id, **details))
        issues, status = [Issue(IType.INFO, f'CREATE acls: The acl created successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE acls: The acl could not be created.')], 500
    return issues, content, status


def read_acls(uuid=None, **kwargs):
    content, count = None, 0
    try:
        content, count = get_query(DBSession, ACL, uuid=uuid, **kwargs)
        # TODO: if chado_id
        if uuid:
            content = content.first()
            content['details'], count = get_query(DBSession, ACLDetail, acl_id=content.id)
            content['details'] = content['details'].all()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ acls: The acls were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ acls: The acls could not be read.')], 500
    return issues, content, count, status


def update_acls(uuid, **kwargs):
    content=None
    try:
        content = get_query(DBSession, ACL, uuid=uuid)[0].update(**kwargs)
        details = kwargs.get('details')
        if isinstance(details, (list, tuple)):
            # Removing missing details
            DBSession.query(ACLDetail).filter(acl_id=content.id).filter(id.out_([d.id for d in details])).delete(synchronize_session='fetch')
            for d in details:
                if 'id' in d:
                    # Updating existing details
                    DBSession.query(ACLDetail).filter(id=d.id).update(**d)
                else:
                    # Adding new details
                    DBSession.add(ACLDetail(acl_id=content.id, **d))
        elif isinstance(details, (dict, ACLDetail)):
            DBSession.add(ACLDetail(acl_id=content.id, **details))
        issues, status = [Issue(IType.INFO, f'UPDATE acls({uuid}): The acls were successfully updated.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'UPDATE acls: The acls could not be updated.')], 500
    return issues, content, status


def delete_acls(uuid=None, **kwargs):
    content = 0
    try:
        content = get_query(DBSession, ACL, uuid=uuid, **kwargs)[0].delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE acls({uuid}): The acls were successfully deleted.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE acls: The acls could not be deleted.')], 500
    return issues, content, status


"""
ObjectTypes service
"""


def read_obj_types(id=None, **kwargs):
    content, count = None, 0
    try:
        content, count = get_query(DBSession, ObjectType, id=id, **kwargs)
        if id or kwargs['value'].get('name'):
            content = content.first()
            content['permission_types'] = get_permission(otype_id=content.id).all()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ object_types: The object_types were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ object_types: The object_types could not be read.')], 500
    return issues, content, count, status


def get_permission(otype_id):
    from biobarcoding.db_models.sysadmin import ObjectTypePermissionType
    return DBSession.query(ObjectTypePermissionType.permission_type).filter(
        ObjectTypePermissionType.object_type_id==otype_id)


"""
PermissionTypes service
"""


def read_perm_types(id=None, uuid=None, type=None, **kwargs):
    content, count = None, 0
    try:
        content, count = get_query(DBSession, PermissionType, id=id, uuid=uuid, name=type, **kwargs)
        if id or type:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ permission_types: The permission_types were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ permission_types: The permission_types could not be read.')], 500
    return issues, content, count, status
