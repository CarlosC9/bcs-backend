"""
A set of hierarchies
Hierarchies have levels
Hierarchies have nodes

Hierarchy=dict(name="", type="", levels=list(dict()), attributes=dict())
HierarchyLevel=dict(name="", hierarchy=<h>, attributes=dict())
HierarchyNode=dict(name="", hierarchy=<h>, level=<l>, parent=<n>, reference_node=<n>, attributes=dict())

"""
from sqlalchemy import and_

from . import get_query
from ..common.pg_helpers import create_or_update_entity
from ..db_models.hierarchies import Hierarchy, HierarchyNode, HierarchyLevel


def list_hierarchies(session, kwargs):
    """
    Enumeration of hierarchies and attributes (id, uuid, name, type -list or hierarchy-)

    :return:
    """
    # TODO Would be good to obtain also the number of levels and the number of nodes in each hierarchy
    query, count = get_query(session, Hierarchy, **kwargs)
    return query.all()


def get_just_hierarchy(session, h_id):
    if not isinstance(h_id, Hierarchy):
        try:
            h = session.query(Hierarchy).get(h_id)
        except:
            h = session.query(Hierarchy).filter(Hierarchy.uuid == h_id).first()
    else:
        h = h_id
    return h


def get_children_of_node(session, node_id):
    _ = []
    children = session.query(HierarchyNode).filter(HierarchyNode.parent_id == node_id).all()
    for child in children:
        d_ = dict(id=child.id, uuid=child.uuid, level_id=child.level_id,
                  reference_node_id=child.reference_node_id, attributes=child.attributes)
        d_["children"] = get_children_of_node(session, child.id)
        _.append(d_)
    return _


def get_hierarchy(session, h_id, format_):
    """
    Obtain the full hierarchy, either flat or nested
    Also the list of levels (without nodes)
    dict(levels=<list of dict(id=<level id>, name=<level name>)>,
         nodes=dict(id=<id>,
                    name="",
                    parent=<id of parent node, if flat>,
                    children=<list of children nodes, if nested>,
                    level_id=<id of hierarchy level>,
                    reference_node_id=<id of referenced node, if it applies>,
                    attributes=<dict>)


    :param session:
    :param h_id:
    :param format_:
    :return:
    """
    h = get_just_hierarchy(session, h_id)
    if h:
        d = dict(id=h.id, uuid=str(h.uuid), type_id=h.h_type_id, name=h.name, attributes=h.attributes, nodes=[], levels=[])
        if format_ == "flat":
            h_nodes = session.query(HierarchyNode).filter(HierarchyNode.hierarchy_id == h.id).all()
            for n in h_nodes:
                d_ = dict(id=n.id, uuid=str(n.uuid), name=n.name, parent_id=n.parent_node_id, level_id=n.level_id,
                          reference_node_id=n.reference_node_id, attributes=n.attributes)
                d["nodes"].append(d_)
        elif format_ == "nested":
            h_nodes = session.query(HierarchyNode).filter(
                and_(HierarchyNode.hierarchy_id == h.id, HierarchyNode.parent_id is None)
            ).all()
            for n in h_nodes:
                d_ = dict(id=n.id, uuid=str(n.uuid), name=n.name, level_id=n.level_id,
                          reference_node_id=n.reference_node_id, attributes=n.attributes)
                d_["children"] = get_children_of_node(session, n.id)
                d["nodes"].append(d_)

        # levels = session.query(HierarchyLevel).filter(HierarchyLevel.hierarchy_id == h.id).all()
        # for level in levels:
        #   level_dict = dict(id=<id>, name=<level name>)
        #   d["levels"].append(level_dict)
    else:
        d = None
    return d


def get_nodes_in_hierarchy_level(session, hl_id):
    """
    Obtain a list of nodes members of hierarchy level "hl_id" (to obtain the levels in a hierarchy, call "get_hierarchy")

    :param session:
    :param hl_id:
    :return:
    """
    pass


# ----------------------------------------------------------------------------------------------------------------------

def create_or_update_hierarchy(session, uuid, name, type_, attributes):
    """
    Create or update a hierarchy. If it exists, do not create again

    :param session:
    :param uuid:
    :param name:
    :param type_:
    :param attributes:
    :return:
    """
    _ = []
    if uuid is not None:
        _.append(("uuid", uuid))
    if name is not None:
        _.append(("name", name))
    if type_ is not None:
        _.append(("h_type_id", type_))
    if attributes is not None:
        _.append(("attributes", attributes))

    _ = zip(*_)
    ent = create_or_update_entity(session, Hierarchy, _[0], _[1], update=True)
    return ent


def create_update_or_delete_hierarchy_levels(session, h_id, levels_list):
    """
    Add or modify levels, delete levels
    Each level is a dict:
      - operation: CU or D
      - name: level name
      - attributes:

    :param h_id:
    :param levels_list:
    :return:
    """


def create_update_or_delete_hierarchy_nodes(session, h_id, node_list):
    """
    Add nodes, delete nodes, modify nodes
    Each node in node_list is a dict:
      - operation: CU or D
      - attributes:

    HierarchyNode=dict(name="", hierarchy=<h>, level=<l>, parent=<n>, reference_node=<n>, attributes=dict())

    :param session:
    :param h_id:
    :param node_list:
    :return:
    """
    objs = dict()  # uuid -> HierarchyNode or Hierarchy or HierarchyLevel
    h = get_just_hierarchy(session, h_id)
    for node in node_list:
        operation = node.get("operation", "cu").lower()
        if operation == "cu":  # Create or update
            _ = [("hierarchy", h)]
            for attr in ["uuid", "name", "level_uuid", "level_id", "parent_uuid", "parent_id",
                         "reference_node_uuid", "reference_node_id", "attributes"]:
                if attr in node:
                    if attr.endswith("_uuid"):
                        v = node[attr]
                        o = objs.get(v)
                        if o is None:
                            if attr == "level_uuid":
                                clazz = HierarchyLevel
                            else:
                                clazz = HierarchyNode
                            o = session.query(clazz).filter(clazz.uuid).first()
                            objs[v] = o
                        _.append((attr[:-len("_uuid")], o))
                    else:
                        _.append((attr, node[attr]))
            _ = zip(*_)
            ent = create_or_update_entity(session, HierarchyNode, _[0], _[1], update=True)
            objs[ent.uuid] = ent
        else:
            # TODO Check if it is a parent node or it is referenced by other nodes.
            #      In that case it should not be deleted
            try:
                o = session.query(HierarchyNode).get(node.get("id"))
            except:
                o = session.query(HierarchyNode).filter(HierarchyNode.uuid == node.get("uuid"))
            o.delete()
