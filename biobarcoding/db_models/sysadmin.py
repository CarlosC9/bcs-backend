import datetime

from sqlalchemy import Column, Integer, String, ForeignKey, DateTime
from sqlalchemy.orm import relationship, backref

from biobarcoding.db_models import ORMBase, GUID, ObjectType


class Format:
    pass


# AUTHENTICATION / AUTHORIZATION

class Authenticator(ORMBase):  # CODES
    """ List of valid authenticators """
    __tablename__ = "authenticators"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    validation_endpoint = Column(String(1024))


class Identity(ORMBase):
    """ Identities.
     Users can have one or more authenticators
     Permissions (ACLs) are assigned to Identities or Groups of identities """
    __versioned__ = {}
    __tablename__ = "identities"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    email = Column(String(80))
    creation_time = Column(DateTime, default=datetime.datetime.utcnow())
    deactivation_time = Column(DateTime)


class IdentityAuthenticator(ORMBase):
    """ Recognized identity authenticator """
    __tablename__ = "identities_authenticators"

    identity_id = Column(Integer, ForeignKey(Identity.id), nullable=False, primary_key=True)
    identity = relationship(Identity, backref=backref("authenticators", cascade="all, delete-orphan"))
    authenticator_id = Column(Integer, ForeignKey(Authenticator.id), nullable=False, primary_key=True)
    authenticator = relationship(Authenticator)

    email = Column(String(80))


class Group(ORMBase):
    __tablename__ = "groups"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class GroupIdentity(ORMBase):
    __tablename__ = "groups_identities"
    group_id = Column(Integer, ForeignKey(Group.id), nullable=False, primary_key=True)
    identity_id = Column(Integer, ForeignKey(Identity.id), nullable=False, primary_key=True)
    group = relationship(Group, backref=backref("identities", cascade="all, delete-orphan"))
    identity = relationship(Identity, backref=backref("groups", cascade="all, delete-orphan"))


class PermissionType(ORMBase):  # CODES
    __tablename__ = "permission_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, nullable=False)  # Object ID (ACL are on objects with UUID)
    name = Column(String(80))


class ACL(ORMBase):
    """ List of permissions on an object. The detail """
    __tablename__ = "permissions"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, nullable=False)  # Object ID (ACL are on objects with UUID)
    object_type = Column(Integer, ForeignKey(ObjectType.id), nullable=False)
    object_id = Column(GUID, nullable=False)


class ACLDetail(ORMBase):
    __tablename__ = "permissions_detail"

    id = Column(Integer, primary_key=True, autoincrement=True)
    acl_id = Column(Integer, ForeignKey(ACL.id))
    acl = relationship(ACL, backref=backref("detail", cascade="all, delete-orphan"))
    identity_id = Column(Integer, ForeignKey(Identity.id))
    group_id = Column(Integer, ForeignKey(Group.id))
    permission_id = Column(Integer, ForeignKey(PermissionType.id))  # Read, Export (Share?), Modify, Delete (depends on the type of object)
    permission = relationship(PermissionType)


# TASKS

class Task:  # Celery task
    """
    Submitted task -To Celery-. Status, log, result, start/end time
    """
    pass
