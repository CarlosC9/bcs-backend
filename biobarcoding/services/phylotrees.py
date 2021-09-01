import os.path
import time

from ..db_models import DBSession as db_session, DBSessionChado as chado_session
from ..db_models.bioinformatics import PhylogeneticTree
from ..db_models.chado import Phylotree, Phylonode
from ..rest import Issue, IType, filter_parse
from . import get_query, get_bioformat


##
# CREATE
##

def create(**kwargs):
    issues = [Issue(IType.WARNING, 'CREATE phylotrees: dummy completed')]
    return issues, None, 200


##
# READ
##

def read(id=None, **kwargs):
    content, count = None, 0
    try:
        content, count = __get_query(id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ phylotrees: The phylotrees were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ phylotrees: The phylotrees could not be read.')], 400
    return issues, content, count, status


##
# UPDATE
##

def update(phylotree_id, **kwargs):
    issues = [Issue(IType.WARNING, 'UPDATE phylotrees: dummy completed')]
    return issues, None, 200


##
# DELETE
##

def __delete_from_bcs(*args):
    db_session.query(PhylogeneticTree).filter(PhylogeneticTree.chado_id.in_(args)) \
        .delete(synchronize_session='fetch')


def delete(id=None, **kwargs):
    content = None
    try:
        content, count = __get_query(id, **kwargs)
        for phylo in content.all():
            __delete_from_bcs(phylo.phylotree_id)
        content = content.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO,
                                f'DELETE phylotrees: The {content} phylotrees were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE phylotrees: The phylotrees could not be removed.')], 404
    return issues, content, status


##
# IMPORT
##

# TODO: NGD newick phylotree import
"""
phylotree:
 - type_id ? cvterm['phylogeny'].cvterms_id
 - for each selected cvterms: new phylotreeprop
 - store the input_file content into a phylotreeprop EDAM_newick
phylonode:
 - type_id ? (root,leaf,internal)
"""
def import_file(input_file, format='newick', **kwargs):
    content = None
    format = get_bioformat(input_file, format or 'newick')
    try:
        # Check phylotree file format
        from Bio import Phylo
        tree = Phylo.read(input_file, format)
    except Exception as e:
        issues = [Issue(IType.ERROR, f'IMPORT phylotress: The file {input_file}.{format} could not be imported.')]
        return issues, content, 409
    try:
        if not kwargs.get('name'):
            kwargs['name'] = os.path.basename(input_file)
        # Create the new phylotree
        phylotree = __new_phylotree(kwargs.get('name'), kwargs.get('comment'), kwargs.get('analysis_id'))
        # Get phylonodes insertion
        phylonodes = __tree2phylonodes(phylotree.phylotree_id, tree.root, None, [0])
        issues, status = [Issue(IType.INFO,
                                'IMPORT phylotrees: The phylotree was successfully imported.',
                                os.path.basename(input_file))], 200
    except Exception as e:
        issues, status = [Issue(IType.ERROR,
                                f'IMPORT phylotress: The file {input_file}.{format} could not be imported.',
                                os.path.basename(input_file))], 409
    return issues, content, status


def __new_phylotree(name, comment=None, analysis_id=None):
    chadotree = __phylotree2chado(name, comment, analysis_id)
    __phylotree2bcs(chadotree)
    return chadotree


def __phylotree2chado(name, comment=None, analysis_id=None):
    dbxref = __new_phylotree_dbxref(name)
    from biobarcoding.db_models.chado import Analysis
    phylotree = Phylotree(dbxref_id=dbxref.dbxref_id, name=name)
    if comment:
        phylotree.comment = comment
    if analysis_id:
        phylotree.analysis_id = chado_session.query(Analysis).filter(
            Analysis.analysis_id == analysis_id).first().analysis_id
    chado_session.add(phylotree)
    chado_session.flush()
    return phylotree


def __new_phylotree_dbxref(name):
    from biobarcoding.db_models.chado import Db, Dbxref
    dbxref = Dbxref(
        db_id=chado_session.query(Db).filter(Db.name == 'null').first().db_id,
        accession=f'phylotree:{name}',
        version=time.strftime("%Y %b %d %H:%M:%S"))
    chado_session.add(dbxref)
    chado_session.flush()
    return dbxref


def __phylotree2bcs(phylotree):
    from biobarcoding.services import get_or_create
    bcs_phylotree = get_or_create(db_session, PhylogeneticTree,
                                  chado_id=phylotree.phylotree_id,
                                  chado_table='phylotree',
                                  name=phylotree.name)
    return bcs_phylotree


def __tree2phylonodes(phylotree_id, node, parent_id=None, index=[0]):
    from biobarcoding.db_models.chado import Phylonode, Feature
    phylonodes = []
    feature_id = chado_session.query(Feature.feature_id).filter(Feature.uniquename == node.name).first()
    phylonode = Phylonode(
        phylotree_id=phylotree_id,
        parent_phylonode_id=parent_id,
        feature_id=feature_id,
        label=node.name,
        distance=node.branch_length,
        left_idx=index[0],
        right_idx=index[0] + 1)
    index[0] += 1
    chado_session.add(phylonode)
    chado_session.flush()
    if len(node.clades) > 0:
        for clade in node.clades:
            phylonodes += __tree2phylonodes(phylotree_id, clade, phylonode.phylonode_id, index)
        phylonode.right_idx = index[0]
        chado_session.add(phylonode)
        chado_session.flush()
    index[0] += 1
    return phylonodes + [phylonode]


##
# EXPORT
##

def export(id=None, format='newick', **kwargs):
    # NGD newick phylotree export
    try:
        if __get_query(id)[0].first():
            __tree2file(id, format, f'/tmp/output_ngd.{format}')
        issues, status = [Issue(IType.INFO, 'EXPORT phylotrees: The phylotree were successfully exported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'EXPORT phylotrees: The phylotree could not be exported.')], 404
    return issues, f'/tmp/output_ngd.{format}', status


def __tree2file(phylotree_id, format, output_file):
    result = ''
    if format == 'newick':
        root = chado_session.query(Phylonode).filter(Phylonode.phylotree_id == phylotree_id,
                                                     Phylonode.parent_phylonode_id == None).first()
        result = __tree2newick(root)
    with open(output_file, "w") as file:
        file.write(result)
    return output_file


def __tree2newick(node, is_root=True):
    result = ''
    children = chado_session.query(Phylonode).filter(Phylonode.parent_phylonode_id == node.phylonode_id)
    if children.count() > 0:
        result += '(' if is_root else '\n('
        i = 1
        for n in children.all():
            result += __tree2newick(n, False)
            result += ',' if i < children.count() else ')'
            i += 1
    if node.label or node.distance:
        label = node.label if node.label else ''
        distance = node.distance if node.distance else '0.00000'
        result += f'\n{label}:{distance}'
    result += ';' if is_root else ''
    return result


##
# GETTER AND OTHERS
##

def __get_query(phylotree_id=None, **kwargs):
    if phylotree_id:
        query = chado_session.query(Phylotree).filter(Phylotree.phylotree_id == phylotree_id)
        return query, query.count()
    phy_clause = {'phylotree_id': {'op': 'in', 'unary': db_session.query(PhylogeneticTree.chado_id).all()}}
    if kwargs.get('filter'):
        try:
            kwargs['filter'] += [phy_clause]
        except Exception as e:
            kwargs['filter'] = [kwargs['filter']] + [phy_clause]
    else:
        kwargs['filter'] = [phy_clause]
    return get_query(chado_session, Phylotree, aux_filter=__aux_own_filter, aux_order=__aux_own_order, **kwargs)


def __aux_own_filter(filter):
    clause = []

    if 'feature_id' in filter:
        _ids = chado_session.query(Phylonode.phylotree_id) \
            .filter(filter_parse(Phylonode, [{'feature_id': filter.get('feature_id')}]))
        clause.append(Phylotree.phylotree_id.in_(_ids))

    if 'organism_id' in filter:
        from biobarcoding.db_models.chado import Feature
        _ids = chado_session.query(Feature.feature_id) \
            .filter(filter_parse(Feature, [{'organism_id': filter.get('organism_id')}]))
        _ids = chado_session.query(Phylonode.phylotree_id) \
            .filter(Phylonode.feature_id.in_(_ids))
        clause.append(Phylotree.phylotree_id.in_(_ids))

    if "prop_cvterm_id" in filter:
        from biobarcoding.db_models.chado import Phylotreeprop
        _ids = chado_session.query(Phylotreeprop.phylotree_id) \
            .filter(filter_parse(Phylotreeprop, [{'type_id': filter.get('prop_cvterm_id')}]))
        clause.append(Phylotree.phylotree_id.in_(_ids))

    from biobarcoding.db_models.chado import Analysis
    if "program" in filter:
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, [{'program': filter.get('program')}]))
        clause.append(Phylotree.analysis_id.in_(_ids))

    if "programversion" in filter:
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, [{'programversion': filter.get('programversion')}]))
        clause.append(Phylotree.analysis_id.in_(_ids))

    if "algorithm" in filter:
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, [{'algorithm': filter.get('algorithm')}]))
        clause.append(Phylotree.analysis_id.in_(_ids))

    from datetime import datetime
    if "added-from" in filter:
        filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, {'timeexecuted': filter.get("added-from")}))
        clause.append(Phylotree.analysis_id.in_(_ids))
    if "added-to" in filter:
        filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, {'timeexecuted': filter.get("added-to")}))
        clause.append(Phylotree.analysis_id.in_(_ids))

    return clause


def __aux_own_order(order):
    # query = query.order(order_parse(Feature, kwargs.get('order'), __aux_own_order))
    return []
