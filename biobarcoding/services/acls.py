from biobarcoding.db_models import DBSession, ObjectType
from biobarcoding.db_models.bioinformatics import BioinformaticObject
from biobarcoding.db_models.sysadmin import ACL, ACLDetail, PermissionType
from biobarcoding.rest import Issue, IType
from biobarcoding.rest import get_query
from biobarcoding.services import get_simple_query, orm2json


def create_acls(**kwargs):
    content = None
    try:
        if kwargs.get('object_uuid'):
            if not kwargs.get('object_type'):
                # TODO: what for the non-bos ?
                kwargs['object_type'] = DBSession.query(BioinformaticObject.bo_type_id)\
                    .filter(BioinformaticObject.uuid==kwargs.get('object_uuid')).first()
        elif kwargs.get('object_type'):
            if not kwargs.get('chado_id'):
                raise
            kwargs['object_uuid'] = DBSession.query(BioinformaticObject)\
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


def read_acls(id=None, uuid=None, **kwargs):
    content, count = None, 0
    try:
        qparams = kwargs['value'] if kwargs.get('value') else kwargs
        # IDs: id, uuid, object_uuid, chado_id + object_type
        if qparams.get('chado_id') and not id and not uuid and not qparams.get('object_uuid'):
            if not qparams.get('object_type'):
                raise
            qparams['object_uuid'] = DBSession.query(BioinformaticObject)\
                .filter(BioinformaticObject.chado_id==qparams.pop('chado_id'),
                        BioinformaticObject.bo_type_id==qparams.get('object_type')).one().uuid
        content, count = get_query(DBSession, ACL, id=id, uuid=uuid, **kwargs)

        if id or uuid or qparams.get('object_uuid'):
            content = orm2json(content.one())
            content['details'] = get_simple_query(DBSession, ACLDetail, acl_id=content.get('id')).all()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ acls: The acls were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ acls: The acls could not be read.')], 500
    return issues, content, count, status


def update_acls(id=None, uuid=None, **kwargs):
    content = None
    try:
        details = kwargs.pop('details')
        content = get_simple_query(DBSession, ACL, id=id, uuid=uuid)
        if kwargs:
            content.update(**kwargs)
        content = content.one()
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
        issues, status = [Issue(IType.INFO, f'UPDATE acls({id}): The acls were successfully updated.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'UPDATE acls: The acls could not be updated.')], 500
    return issues, content, status


def delete_acls(id=None, uuid=None, **kwargs):
    content = 0
    try:
        content = get_query(DBSession, ACL, id=id, uuid=uuid, **kwargs)[0].delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE acls({id}): The acls were successfully deleted.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE acls: The acls could not be deleted.')], 500
    return issues, content, status


"""
ObjectTypes service
"""


def read_obj_types(id=None, uuid=None, **kwargs):
    content, count = None, 0
    try:
        content, count = get_query(DBSession, ObjectType, id=id, uuid=uuid, **kwargs)
        if id or uuid:
            content = orm2json(content.first())
            from biobarcoding.db_models.sysadmin import ObjectTypePermissionType
            perms = DBSession.query(PermissionType).join(ObjectTypePermissionType)\
                .filter(ObjectTypePermissionType.object_type_id==content.get('id'))
            content['permission_types'] = perms.all()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ object_types: The object_types were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ object_types: The object_types could not be read.')], 500
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
        issues, status = [Issue(IType.ERROR, 'READ permission_types: The permission_types could not be read.')], 500
    return issues, content, count, status
