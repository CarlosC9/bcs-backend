"""
Functions to determine if ACL objects can be accessed for certain operation(s)

* check if object can be acted by user to perform an operation -> True/False
* check which operations can perform a user on an object -> List of operations
* prepare WHERE subclause to speed up looking up for ACL objects, given user and an operation (or list of operations) -> SQL Where Clause

"""
from typing import Set
from datetime import datetime

from sqlalchemy import and_, or_

from ..db_models.bioinformatics import BioinformaticObject
from ..db_models.sysadmin import Identity, ACLExpression, ACL, PermissionType, ACLDetail


def check(sess, o, ident: Identity, permission: PermissionType) -> bool:
    """

    :param sess: Database session
    :param o: Object to check
    :param ident:
    :param permission:
    :return:
    """

    if isinstance(o, BioinformaticObject):
        obj_type_id = o.bo_type_id
        o_uuid = o.uuid

    ahora = datetime.now()
    # TODO Fast method
    #  - Find groups, organizations and roles to which the current user is member of: {groups(ident)}, {roles(ident)}, {organizations(ident)}
    #  - an object can have ACLdetails directly assigned or inherited from: object type, collection(s) to which it is member, others
    #  - Search: valid ACLDetails assigned to "o_uuid", for "permission",
    #            group IN {groups(ident)} OR
    #            role IN {roles(ident)} OR
    #            organizations IN {organizations(ident)}
    #  Permission granted if any record found. Positive permission
    acl_details = sess.query(ACLDetail). \
        filter(and_(
        or_(ACLDetail.validity_start is None, ACLDetail.validity_start <= ahora),
        or_(ACLDetail.validity_end is None, ACLDetail.validity_end > ahora))). \
        join(ACLDetail.acl).filter(
        and_(ACL.object_type == obj_type_id, ACL.object_uuid == o_uuid)).first()

    # Rule stored in the database. Find the active one for the desired function
    # TODO if there is no ACLExpression, search ACL elements (or maybe generate an expression from the ACL elements)
    rule = sess.query(ACLExpression.expression). \
        filter(and_(
        or_(ACLExpression.validity_start is None, ACLExpression.validity_start <= ahora),
        or_(ACLExpression.validity_end is None, ACLExpression.validity_end > ahora))). \
        join(ACLExpression.acl).filter(
        and_(ACL.object_type == obj_type_id, ACL.object_uuid == o_uuid)).first()


def check_operations(o, user) -> Set[str]:
    pass


def get_where_subclause(user) -> str:
    # TODO Fast method
    #  - Find groups, organizations and roles to which the current user is member of: {groups(ident)}, {roles(ident)}, {organizations(ident)}
    #  - an object can have ACLdetails directly assigned or inherited from: object type, collection(s) to which it is member, others
    #  - Search: valid ACLDetails assigned to "o_uuid", for "permission",
    #            group IN {groups(ident)} OR
    #            role IN {roles(ident)} OR
    #            organizations IN {organizations(ident)}
    #  Permission granted if any record found. Positive permission

    pass
