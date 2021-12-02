import uuid as uuid
from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship, backref

from . import ORMBase, GUID

prefix = "hie_"


class HierarchyType(ORMBase):
    __tablename__ = f"{prefix}hierarchy_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))


class Hierarchy(ORMBase):
    __tablename__ = f"{prefix}hierarchies"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    h_type_id = Column(Integer, ForeignKey(HierarchyType.id))
    h_type = relationship(HierarchyType)
    attributes = Column(JSONB)


class HierarchyLevel(ORMBase):
    __tablename__ = f"{prefix}hierarchy_levels"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(128))
    uuid = Column(GUID, unique=True, default=uuid.uuid4)

    hierarchy_id = Column(Integer, ForeignKey(Hierarchy.id), nullable=False, primary_key=False)
    hierarchy = relationship(Hierarchy, backref=backref("levels", cascade="all, delete-orphan"))
    attributes = Column(JSONB)


class HierarchyNode(ORMBase):
    __tablename__ = f"{prefix}hierarchy_nodes"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(128))
    uuid = Column(GUID, unique=True, default=uuid.uuid4)

    hierarchy_id = Column(Integer, ForeignKey(Hierarchy.id), nullable=False, primary_key=False)
    # Backref -backref=backref("nodes", cascade="all, delete-orphan")- not used because a hierarchy can contain millions of nodes
    hierarchy = relationship(Hierarchy)

    parent_node_id = Column(Integer, ForeignKey(f"{prefix}hierarchy_nodes.id"), nullable=True, primary_key=False)
    # Backref -backref=backref("children", cascade="all, delete-orphan")- not used because a node may be thousands of children?
    parent_node = relationship("HierarchyNode", foreign_keys=[parent_node_id], remote_side=[id])

    level_id = Column(Integer, ForeignKey(HierarchyLevel.id), nullable=True, primary_key=False)
    # Backref -backref=backref("nodes", cascade="all, delete-orphan")- not used because a hierarchy can contain millions of nodes
    level = relationship(HierarchyLevel)

    reference_node_id = Column(Integer, ForeignKey(f"{prefix}hierarchy_nodes.id"), nullable=True, primary_key=False)
    reference_node = relationship("HierarchyNode", foreign_keys=[reference_node_id], remote_side=[id])

    attributes = Column(JSONB)
