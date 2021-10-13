import os.path
import time

from .ontologies import get_cvterm_query
from ..db_models import DBSession as db_session, DBSessionChado as chado_session
from ..db_models.bioinformatics import PhylogeneticTree, bio_object_type_id
from ..db_models.chado import Phylotree, Phylonode, Feature
from ..rest import Issue, IType, filter_parse, auth_filter
from . import get_query, get_bioformat, get_or_create, get_orm_params


##
# CREATE
##

def __check_phylo_values(**values):
    if not values.get('dbxref_id'):
        from biobarcoding.db_models.chado import Db, Dbxref
        values['dbxref_id'] = get_or_create(chado_session, Dbxref,
                                            db_id=chado_session.query(Db).filter(Db.name == 'null').one().db_id,
                                            accession=f'phylotree:{values.get("name")}',
                                            version=time.strftime("%Y %b %d %H:%M:%S")).dbxref_id

    if not values.get('type_id') and values.get('type'):
        try:
            values['type_id'] = get_cvterm_query(type=values.get('type'), subtype=values.get('subtype'))[0].one().cvterm_id
        except Exception as e:
            pass
    if not values.get('type_id'):
        try:
            values['type_id'] = get_cvterm_query(type='phylotree')[0].one().cvterm_id
        except Exception as e:
            pass

    # Analysis row could exist for jobs, so get or create
    ansis_trigger = ('job_id', 'program', 'programversion', 'sourcename')
    if not values.get('analysis_id') and any(key in values.keys() for key in ansis_trigger):
        from .analyses import __get_query as get_ansis_query, create as create_ansis
        try:
            unique_keys = ['job_id'] if values.get('job_id') else ('program', 'programversion', 'sourcename')
            content, count = get_ansis_query(**{k:values[k] for k in unique_keys if k in values})
            content = content.one()
        except Exception as e:
            treename, values['name'] = values.get('name'), None     # Avoid to spread the treename to analysis
            issues, content, status = create_ansis(**values)
            values['name'] = treename
            chado_session.flush()   # analysis_id required
        values['analysis_id'] = content.analysis_id

    return get_orm_params(Phylotree, **values)


def create(**kwargs):
    content = None
    try:
        values = __check_phylo_values(**kwargs)
        content = Phylotree(**values)
        chado_session.add(content)
        # tree to bcs
        chado_session.flush()   # phylotree_id required
        phylo = get_or_create(db_session, PhylogeneticTree,
                              chado_id=content.phylotree_id,
                              chado_table='phylotree',
                              name=content.name)
        issues, status = [Issue(IType.INFO,
                                f'CREATE phylotrees: The phylotree "{kwargs.get("name")}" created successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR,
                                f'CREATE phylotrees: The phylotree "{kwargs.get("name")}" could not be created.')], 409
    return issues, content, status


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

def __delete_from_bcs(*ids):
    query = db_session.query(PhylogeneticTree).filter(PhylogeneticTree.chado_id.in_(ids))
    return query.count()
    # TODO: check why all rows are deleted in bcs without filtering
    # return query.delete(synchronize_session='fetch')


def delete(id=None, **kwargs):
    content = None
    try:
        content, count = __get_query(id, purpose='delete', **kwargs)
        ids = [phylo.phylotree_id for phylo in content.all()]
        bcs_delete = __delete_from_bcs(*ids)
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

def import_file(input_file, format=None, **kwargs):
    ## TODO: NGD newick phylotree import
    # phylotree: store the input_file content into a phylotreeprop EDAM_newick ?
    # phylonode: type_id ? (root,leaf,internal)
    content = None
    format = get_bioformat(input_file, format)
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
        issues, content, status = create(**kwargs)
        # Get phylonodes insertion
        phylonodes = __tree2phylonodes(content.phylotree_id, tree.root, None, [0])
        issues, status = [Issue(IType.INFO,
                                'IMPORT phylotrees: The phylotree was successfully imported.',
                                os.path.basename(input_file))], 200
    except Exception as e:
        issues, status = [Issue(IType.ERROR,
                                f'IMPORT phylotress: The file {input_file}.{format} could not be imported.',
                                os.path.basename(input_file))], 409
    return issues, content, status


def __tree2phylonodes(phylotree_id, node, parent_id=None, index=[0]):
    phylonodes = []
    feature_id = chado_session.query(Feature.feature_id).filter(Feature.uniquename == node.name).first()
    phylonode = get_or_create(chado_session, Phylonode,
                              phylotree_id=phylotree_id,
                              parent_phylonode_id=parent_id,
                              feature_id=feature_id,
                              label=node.name,
                              distance=node.branch_length,
                              left_idx=index[0],
                              right_idx=index[0] + 1)
    # Check for children
    index[0] += 1
    if len(node.clades) > 0:
        for clade in node.clades:
            phylonodes += __tree2phylonodes(phylotree_id, clade, phylonode.phylonode_id, index)
        phylonode.right_idx = index[0]
        chado_session.merge(phylonode)
    index[0] += 1
    return [phylonode] + phylonodes


##
# EXPORT
##

def export(id=None, format='newick', **kwargs):
    # NGD newick phylotree export
    try:
        if id and __get_query(id, purpose='share')[0].one():
            __tree2file(id, format, f'/tmp/output_ngd.{format}')
        else:
            # TODO: allow export by filter ?
            raise Exception('EXPORT phylotrees: The phylotree is missing to export.')
        issues, status = [Issue(IType.INFO, 'EXPORT phylotrees: The phylotree were successfully exported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'EXPORT phylotrees: The phylotree could not be exported.')], 404
    return issues, f'/tmp/output_ngd.{format}', status


def __tree2file(phylotree_id, format, output_file):
    result = ''
    format = get_bioformat(output_file, format)
    # TODO: build Bio.Phylo.BaseTree from chado, and Phylo.write(tree, output_file, format)
    if format == 'newick':
        root = chado_session.query(Phylonode).filter(Phylonode.phylotree_id == phylotree_id,
                                                     Phylonode.parent_phylonode_id == None).one()
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

def __get_query(phylotree_id=None, purpose='read', **kwargs):
    if phylotree_id:
        query = chado_session.query(Phylotree).filter(Phylotree.phylotree_id == phylotree_id)
        return query, query.count()
    from biobarcoding.db_models.sysadmin import PermissionType
    purpose_id = db_session.query(PermissionType.id).filter(PermissionType.name==purpose).one()
    phy_clause = db_session.query(PhylogeneticTree.chado_id) \
        .filter(auth_filter(PhylogeneticTree, purpose_id, [bio_object_type_id['phylogenetic-tree']]))
    phy_clause = [i for i, in phy_clause.all()]
    phy_clause = Phylotree.phylotree_id.in_(phy_clause)
    query = chado_session.query(Phylotree).filter(phy_clause)
    return get_query(chado_session, Phylotree, query=query, aux_filter=__aux_own_filter, aux_order=__aux_own_order, **kwargs)


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
