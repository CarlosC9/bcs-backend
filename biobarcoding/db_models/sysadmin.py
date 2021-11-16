import datetime
import uuid

from sqlalchemy import Column, ForeignKey, UniqueConstraint, Boolean, Integer, BigInteger, String, DateTime, Text, JSON, \
    Sequence
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship, backref

from . import ORMBase, GUID, ObjectType
from .bioinformatics import BioinformaticObject

# AUTHENTICATION / AUTHORIZATION

prefix = "sa_auth_"

authorizable_type_id = {
    "identity": 1,
    "organization": 2,
    "group": 3,
    "role": 4,
}


class Authenticator(ORMBase):  # CODES
    """ List of valid authenticators """
    __tablename__ = f"{prefix}authenticators"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    validation_endpoint = Column(String(1024))


class Authorizable(ORMBase):
    """ Authorizables.
     Permissions (ACLs) are assigned to Authorizables """
    __versioned__ = {}
    __tablename__ = f"{prefix}authorizables"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    authorizable_type_id = Column(Integer, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity': 'authorizable',
        'polymorphic_on': authorizable_type_id
    }


class Identity(Authorizable):
    """ Identities.
     Users can have one or more authenticators
     Permissions (ACLs) are assigned to Identities or Groups of identities """
    __versioned__ = {}
    __tablename__ = f"{prefix}identities"
    __mapper_args__ = {
        'polymorphic_identity': authorizable_type_id['identity'],
    }

    id = Column(Integer, ForeignKey(Authorizable.id), primary_key=True)
    name = Column(String(255))
    email = Column(String(255))
    # configuration = Column(JSON)  # To store personal preferences, from language to other visualization features
    creation_time = Column(DateTime, default=datetime.datetime.utcnow())
    deactivation_time = Column(DateTime)


class IdentityAuthenticator(ORMBase):
    """ Recognized identity authenticator """
    __tablename__ = f"{prefix}identities_authenticators"

    identity_id = Column(Integer, ForeignKey(Identity.id), nullable=False, primary_key=True)
    identity = relationship(Identity, backref=backref("authenticators", cascade="all, delete-orphan"))
    authenticator_id = Column(Integer, ForeignKey(Authenticator.id), nullable=False, primary_key=True)
    authenticator = relationship(Authenticator)

    email = Column(String(255))
    name = Column(String(255))
    authenticator_info = Column(JSON)


class Organization(Authorizable):
    __tablename__ = f"{prefix}organizations"
    __mapper_args__ = {
        'polymorphic_identity': authorizable_type_id['organization'],
    }

    id = Column(Integer, ForeignKey(Authorizable.id), primary_key=True)
    name = Column(String(80))


class OrganizationIdentity(ORMBase):
    __tablename__ = f"{prefix}organizations_identities"
    organization_id = Column(Integer, ForeignKey(Organization.id), nullable=False, primary_key=True)
    identity_id = Column(Integer, ForeignKey(Identity.id), nullable=False, primary_key=True)
    organization = relationship(Organization, backref=backref("identities", cascade="all, delete-orphan"))
    identity = relationship(Identity, backref=backref("organizations", cascade="all, delete-orphan"))


class Group(Authorizable):
    __tablename__ = f"{prefix}groups"
    __mapper_args__ = {
        'polymorphic_identity': authorizable_type_id['group'],
    }

    id = Column(Integer, ForeignKey(Authorizable.id), primary_key=True)
    name = Column(String(80))


class GroupIdentity(ORMBase):
    __tablename__ = f"{prefix}groups_identities"
    group_id = Column(Integer, ForeignKey(Group.id), nullable=False, primary_key=True)
    identity_id = Column(Integer, ForeignKey(Identity.id), nullable=False, primary_key=True)
    group = relationship(Group, backref=backref("identities", cascade="all, delete-orphan"))
    identity = relationship(Identity, backref=backref("groups", cascade="all, delete-orphan"))


class Role(Authorizable):
    __tablename__ = f"{prefix}roles"
    __mapper_args__ = {
        'polymorphic_identity': authorizable_type_id['role'],
    }

    id = Column(Integer, ForeignKey(Authorizable.id), primary_key=True)
    name = Column(String(80))


class RoleIdentity(ORMBase):
    __tablename__ = f"{prefix}roles_identities"
    role_id = Column(Integer, ForeignKey(Role.id), nullable=False, primary_key=True)
    identity_id = Column(Integer, ForeignKey(Identity.id), nullable=False, primary_key=True)
    role = relationship(Role, backref=backref("identities", cascade="all, delete-orphan"))
    identity = relationship(Identity, backref=backref("roles", cascade="all, delete-orphan"))


class SystemFunction(ORMBase):
    """
    Functions of the system, for two purposes:
    * Backend functions: to control execution permissions
    * Frontend functions: to control display permissions (to dynamically prepare the view)

    Backend functions can be directly annotated, they do not need to appear here

    """
    __tablename__ = f"{prefix}functions"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(512))


class PermissionType(ORMBase):  # CODES
    __tablename__ = f"{prefix}permission_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, nullable=False, default=uuid.uuid4)  # Object ID (ACL are on objects with UUID)
    name = Column(String(80))
    rank = Column(Integer, nullable=False, default=0)


class ObjectTypePermissionType(ORMBase):
    __tablename__ = f"{prefix}obj_types_perm_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    object_type_id = Column(Integer, ForeignKey(ObjectType.id), nullable=False, primary_key=True)
    permission_type_id = Column(Integer, ForeignKey(PermissionType.id), nullable=False, primary_key=True)
    object_type = relationship(ObjectType, backref=backref("permission_types", cascade="all, delete-orphan"))
    permission_type = relationship(PermissionType, backref=backref("object_types", cascade="all, delete-orphan"))


class Collection(ORMBase):
    """ List of objects """
    __tablename__ = f"{prefix}collections"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, nullable=False, default=uuid.uuid4)  # Object ID (Collection are on objects with UUID)
    name = Column(String(512))


class CollectionDetail(ORMBase):
    """ Association of an object with a collection """
    __tablename__ = f"{prefix}collections_detail"

    id = Column(Integer, primary_key=True, autoincrement=True)
    collection_id = Column(Integer, ForeignKey(Collection.id), nullable=False)
    collection = relationship(Collection, backref=backref("detail", cascade="all, delete-orphan"))
    object_type = Column(Integer, ForeignKey(ObjectType.id), nullable=False)
    object_uuid = Column(GUID, nullable=False)


class ACL(ORMBase):
    """ List of permissions on an object. The detail """
    __tablename__ = f"{prefix}permissions"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, nullable=False, default=uuid.uuid4)  # Object ID (ACL are on objects with UUID)
    object_type = Column(Integer, ForeignKey(ObjectType.id), nullable=False)
    object_uuid = Column(GUID, nullable=False)

    __table_args__ = (
        UniqueConstraint(object_uuid, object_type, name=__tablename__ + '_c1'),
    )


class ACLExpression(ORMBase):
    """ Authorization in the form of expression, evaluated by """
    __tablename__ = f"{prefix}permissions_expression"

    id = Column(Integer, primary_key=True, autoincrement=True)
    acl_id = Column(Integer, ForeignKey(ACL.id))
    acl = relationship(ACL, backref=backref("expression", cascade="all, delete-orphan"))

    expression = Column(String(500))

    validity_start = Column(DateTime, nullable=True)
    validity_end = Column(DateTime, nullable=True)


class ACLDetail(ORMBase):
    """ Detail can be directly specified or generated by an ACLExpression (which has to be compiled) """
    __tablename__ = f"{prefix}permissions_detail"

    id = Column(Integer, primary_key=True, autoincrement=True)
    acl_id = Column(Integer, ForeignKey(ACL.id), nullable=False)
    acl = relationship(ACL, backref=backref("details", cascade="all, delete-orphan"))

    # Origin of this record (can be NULL if directly specified)
    acl_expression_id = Column(Integer, ForeignKey(ACLExpression.id), nullable=True)
    acl_expression = relationship(ACLExpression, backref=backref("compiled", cascade="all, delete-orphan"))

    authorizable_id = Column(Integer, ForeignKey(Authorizable.id, ondelete="CASCADE"), nullable=False)
    authorizable = relationship(Authorizable)

    permission_id = Column(Integer, ForeignKey(PermissionType.id), nullable=False)  # Read, Export (Share?), Modify, Delete (depends on the type of object)
    permission = relationship(PermissionType)

    validity_start = Column(DateTime, nullable=True)
    validity_end = Column(DateTime, nullable=True)


# TASKS

prefix = "sa_task_"


class TaskStatus(ORMBase):
    __tablename__ = f"{prefix}statuses"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, nullable=False, default=uuid.uuid4)  # Object ID
    name = Column(String(80))


class Task(ORMBase):  # Celery task
    """
    Submitted task -To Celery-. Status, log, result, start/end time
    """
    __tablename__ = f"{prefix}instances"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, nullable=False, default=uuid.uuid4)  # Object ID (ACL are on objects with UUID)
    description = Column(String(80))
    params = Column(JSON)
    creation_time = Column(DateTime, default=datetime.datetime.utcnow())
    finalization_time = Column(DateTime)
    status = Column(Integer, ForeignKey(TaskStatus.id))
    log = Column(Text)


# BROWSER FILTER

prefix = "sa_browser_"


# class BrowserFilterForm(ORMBase):
#     __tablename__ = f"{prefix}filter_forms"
#
#     id = Column(Integer, ForeignKey(ObjectType.id), primary_key=True)
#     # id = Column(Integer, primary_key=True, autoincrement=True)
#     # bo_type_id = Column(Integer, ForeignKey(ObjectType.id), nullable=False, unique=True)
#     uuid = Column(GUID, unique=True)
#     form = Column(Text)


class BrowserFilter(ORMBase):
    __tablename__ = f"{prefix}filters"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80), nullable=False)
    type = Column(String(80), nullable=False)
    values = Column(JSON)
    user_id = Column(Integer, ForeignKey(Identity.id, ondelete="CASCADE"), nullable=False)

    __table_args__ = (
        UniqueConstraint(name, type, user_id, name=__tablename__+'_c1'),
    )


# ANNOTATION FORM

prefix = "sa_annotation_"


class AnnotationFormItem(ORMBase):
    __tablename__ = f"{prefix}form_item"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    name = Column(String(80), nullable=False, unique=True)
    description = Column(String(512))
    cvterm_id = Column(Integer, unique=True)
    dbxref_id = Column(Integer)


class AnnotationFormItemObjectType(ORMBase):
    __tablename__ = f"{prefix}form_item_object_type"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    form_item_id = Column(Integer, ForeignKey(AnnotationFormItem.id), nullable=False)
    object_type_id = Column(Integer, ForeignKey(ObjectType.id), nullable=False)
    form_item = relationship(AnnotationFormItem, backref=backref("annotation_form_items", cascade="all, delete-orphan"))
    object_type = relationship(ObjectType, backref=backref("object_types", cascade="all, delete-orphan"))

    __table_args__ = (
        UniqueConstraint(form_item_id, object_type_id, name=__tablename__ + '_c1'),
    )


class AnnotationFormField(AnnotationFormItem):
    __tablename__ = f"{prefix}form_field"

    id = Column(BigInteger, ForeignKey(AnnotationFormItem.id, ondelete="CASCADE"), primary_key=True)
    type = Column(String(32), nullable=False)   # tag, attribute, relationship
    range = Column(JSONB)
    view_type = Column(String(32))  # check, radio, date, etc.
    multiple = Column(Boolean, default=False)

    __mapper_args__ = {
        'polymorphic_identity': 'annotation_form_field',
        'polymorphic_on': type
    }


class AnnotationFormTag(AnnotationFormField):
    __tablename__ = f"{prefix}form_tag"
    __mapper_args__ = {
        'polymorphic_identity': 'tag',
    }
    id = Column(BigInteger, ForeignKey(AnnotationFormField.id, ondelete="CASCADE"), primary_key=True)


class AnnotationFormAttribute(AnnotationFormField):
    __tablename__ = f"{prefix}form_attribute"
    __mapper_args__ = {
        'polymorphic_identity': 'attribute',
    }
    id = Column(BigInteger, ForeignKey(AnnotationFormField.id, ondelete="CASCADE"), primary_key=True)


class AnnotationFormRelationship(AnnotationFormField):
    __tablename__ = f"{prefix}form_relationship"
    __mapper_args__ = {
        'polymorphic_identity': 'relationship',
    }
    id = Column(BigInteger, ForeignKey(AnnotationFormField.id, ondelete="CASCADE"), primary_key=True)


class AnnotationFormTemplate(AnnotationFormItem):
    __tablename__ = f"{prefix}form_template"

    id = Column(BigInteger, ForeignKey(AnnotationFormItem.id, ondelete="CASCADE"), primary_key=True)


class AnnotationFormTemplateField(ORMBase):
    __tablename__ = f"{prefix}form_template_field"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    form_field_id = Column(Integer, ForeignKey(AnnotationFormField.id), nullable=False)
    form_template_id = Column(Integer, ForeignKey(AnnotationFormTemplate.id), nullable=False)
    form_field = relationship(AnnotationFormField, backref=backref("annotation_form_fields", cascade="all, delete-orphan"))
    form_template = relationship(AnnotationFormTemplate, backref=backref("annotation_form_templates", cascade="all, delete-orphan"))
    name = Column(String(80))
    rank = Column(Integer, Sequence('annotation_form_rank_seq'), nullable=False)
    # rank = Column(Integer, nullable=False, default=select([func.max(1, func.max(rank))]))

    __table_args__ = (
        UniqueConstraint(form_template_id, rank, name=__tablename__ + '_c1'),
    )


# ANNOTATION INSTANCE

class AnnotationItem(ORMBase):
    __tablename__ = f"{prefix}item"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    object_uuid = Column(GUID, ForeignKey(BioinformaticObject.uuid, ondelete="CASCADE"), nullable=False)
    type = Column(String(80), nullable=False)    # text, form, field
    name = Column(String(80), nullable=False)
    rank = Column(Integer, Sequence('annotation_rank_seq'), nullable=False)
    # rank = Column(Integer, nullable=False, default=select([func.max(1, func.max(rank))]))

    __table_args__ = (
        UniqueConstraint(object_uuid, rank, name=__tablename__ + '_c1'),
        UniqueConstraint(object_uuid, name, name=__tablename__ + '_c2'),
    )
    __mapper_args__ = {
        'polymorphic_identity': 'annotation_item',
        'polymorphic_on': type
    }


class AnnotationText(AnnotationItem):
    __tablename__ = f"{prefix}text"
    __mapper_args__ = {
        'polymorphic_identity': 'text',
    }
    id = Column(BigInteger, ForeignKey(AnnotationItem.id, ondelete="CASCADE"), primary_key=True)
    value = Column(String(512))


class AnnotationTemplate(AnnotationItem):
    __tablename__ = f"{prefix}template"
    __mapper_args__ = {
        'polymorphic_identity': 'template',
    }
    id = Column(BigInteger, ForeignKey(AnnotationItem.id, ondelete="CASCADE"), primary_key=True)
    form_template_id = Column(Integer, ForeignKey(AnnotationFormTemplate.id), nullable=False)
    form_template = relationship(AnnotationFormTemplate, backref=backref("annotation_templates", cascade="all, delete-orphan"))
    value = Column(JSONB)


class AnnotationField(AnnotationItem):
    __tablename__ = f"{prefix}field"
    __mapper_args__ = {
        'polymorphic_identity': 'field',
    }
    id = Column(BigInteger, ForeignKey(AnnotationItem.id, ondelete="CASCADE"), primary_key=True)
    form_field_id = Column(Integer, ForeignKey(AnnotationFormField.id), nullable=False)
    form_field = relationship(AnnotationFormField, backref=backref("annotation_fields", cascade="all, delete-orphan"))
    value = Column(JSONB)
