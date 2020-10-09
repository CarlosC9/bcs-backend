from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.orm import relationship, backref

from biobarcoding.db_models import ORMBase, GUID


prefix = "hie_"


class HierarchyType(ORMBase):
    __tablename__ = f"{prefix}_hierarchy_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class Hierarchy(ORMBase):
    __tablename__ = f"{prefix}_hierarchies"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    h_type = Column(Integer, ForeignKey(HierarchyType.id))


class HierarchyLevel(ORMBase):
    __tablename__ = f"{prefix}_hierarchy_levels"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(128))

    hierarchy_id = Column(Integer, ForeignKey(Hierarchy.id), nullable=False, primary_key=False)
    hierarchy = relationship(Hierarchy, backref=backref("levels", cascade="all, delete-orphan"))


class HierarchyNode(ORMBase):
    __tablename__ = f"{prefix}_hierarchy_levels"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(128))

    hierarchy_id = Column(Integer, ForeignKey(Hierarchy.id), nullable=False, primary_key=False)

    parent_node_id = Column(Integer, ForeignKey("HierarchyNode.id"), nullable=False, primary_key=False)
    # Backref -backref=backref("children", cascade="all, delete-orphan")- not used because a node may be thousands of children?
    parent_node = relationship("HierarchyNode")

    # Backref -backref=backref("nodes", cascade="all, delete-orphan")- not used because a hierarchy can contains millions of nodes
    # hierarchy = relationship(Hierarchy)
    level_id = Column(Integer, ForeignKey(HierarchyLevel.id), nullable=False, primary_key=False)
    # Backref -backref=backref("nodes", cascade="all, delete-orphan")- not used because a hierarchy can contains millions of nodes
    level = relationship(HierarchyLevel)

    referred_node_id = Column(Integer, ForeignKey("HierarchyNode.id"), nullable=False, primary_key=False)
    referred_node = relationship("HierarchyNode")

