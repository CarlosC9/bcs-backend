import dendropy
from dendropy.datamodel.treemodel import Tree
from dendropy.datamodel.treecollectionmodel import TreeList
import sys


def get_posterior_height_tree(tree: Tree):
    for n in tree.preorder_node_iter():
        if not n.is_leaf():
            posterior = n.annotations.get_value(name="posterior")
            if posterior:
                posterior = '{:.2f}'.format(float(posterior))
            height = n.annotations.get_value(name="height")
            if height:
                height = '{:.2f}'.format(float(height) / 1e6)

            n.label = f"{height} / {posterior}"

    return tree.as_string("newick")


def get_height_tree(tree: dendropy.datamodel.treemodel.Tree):
    for n in tree.preorder_node_iter():
        if not n.is_leaf():
            height = n.annotations.get_value(name="height")
            if height:
                height = '{:.2f}'.format(float(height) / 1e6)
                n.label = str(height)

    return tree.as_string("newick")


def get_hpd_height_tree(tree: dendropy.datamodel.treemodel.Tree):
    for n in tree.preorder_node_iter():
        if not n.is_leaf():
            hpd = n.annotations.get_value(name="height_95%_HPD")
            if hpd:
                hpd1 = '{:.2f}'.format(float(hpd[0]))
                hpd2 = '{:.2f}'.format(float(hpd[1]))
                n.label = f"{hpd1} - {hpd2}"

    return tree.as_string("newick")


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    ds = dendropy.DataSet.get_from_path(input_filename, "nexus", extract_comment_metadata=True)
    tree = ds.tree_lists[0][0]

    tree_list = TreeList()
    tree_list.read(data=get_posterior_height_tree(tree), schema="newick")
    tree_list[-1].label = "height / posterior"
    tree_list.read(data=get_height_tree(tree), schema="newick")
    tree_list[-1].label = "height"
    tree_list.read(data=get_hpd_height_tree(tree), schema="newick")
    tree_list[-1].label = "height 95% HPD"

    tree_list.write(path=output_filename, schema="nexus")