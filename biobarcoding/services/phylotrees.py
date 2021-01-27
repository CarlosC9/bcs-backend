import os, time

from biobarcoding.authentication import bcs_session
from biobarcoding.db_models import DBSession as db_session
from biobarcoding.db_models import DBSessionChado as chado_session


@bcs_session(read_only=False)
def create_phylotrees(name = None, comment = None):
    return {'status':'success','message':'CREATE: phylotrees dummy completed'}, 200


@bcs_session(read_only=True)
def read_phylotrees(id = None, analysis_id = None, name = None, comment = None, feature_id = None):
    query = __get_query(id, analysis_id, name, comment, feature_id)
    try:
        from biobarcoding.services import chado2json
        if id:
            return chado2json(query)[0], 200
        return chado2json(query), 200
    except Exception as e:
        return f'Unable to get any result for the query.', 500


@bcs_session(read_only=False)
def update_phylotrees(phylotree_id, name = None, comment = None, analysis_id = None):
    return {'status':'success','message':'UPDATE: phylotrees dummy completed'}, 200


def __delete_from_bcs(ids):
    from biobarcoding.db_models.bioinformatics import PhylogeneticTree
    db_session.query(PhylogeneticTree).filter(PhylogeneticTree.chado_phylotree_id.in_(ids))\
        .delete(synchronize_session='fetch')


@bcs_session(read_only=False)
def delete_phylotrees(id=None, ids=None):
    from biobarcoding.db_models.chado import Phylotree, Dbxref
    result = chado_session.query(Phylotree)
    msg='[ '
    if id:
        result = result.filter(Phylotree.phylotree_id==id)
        __delete_from_bcs([id])
        msg+=f'{id} '
    elif ids:
        result = result.filter(Phylotree.phylotree_id.in_(ids))
        __delete_from_bcs(ids)
        msg+=f'{ids} '
    else:
        return {'status': 'failure', 'message': f'Phylotrees IDs missed.'}, 400
    msg+=']'
    result.delete(synchronize_session='fetch')
    # try:
    #     chado_session.commit()
    #     db_session.commit()
    # except Exception as e:
    #     chado_session.rollback()
    #     db_session.rollback()
    #     return {'status':'failure','message':f'Phylotrees could not be deleted.'}, 400
    return {'status':'success','message':f'Phylotrees deleted: {msg}'}, 200


# Legacy import with python-chado
def __import_phylotrees(input_file, name = None, comment = None, analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    try:
        if not analysis_id:
            analysis_id = conn.analysis.get_analyses()[0]['analysis_id']
        response = conn.phylogeny.load_tree(input_file, analysis_id, name=name)
        return {'status':'success','message':f'{response} phylotrees were successfully imported.'}, 200
    except Exception as e:
        print(e)
        return {'status':'failure','message':f'The phylotree could not be imported.'}, 500


@bcs_session(read_only=True)
def export_phylotrees(phylotree_id = None):
    filepath = '/tmp/output_seqs.fas'
    query = __get_query(phylotree_id)
    with open(filepath, "w") as file:
        from biobarcoding.services import chado2json
        file.write(f'{chado2json(query)}')
    return filepath, 200


# NGD newick phylotree import
# TODO:
"""
phylotree:
 - type_id ? cvterm['phylogeny'].cvterms_id
 - for each selected cvterms: new phylotreeprop
 - store the input_file content into a phylotreeprop EDAM_newick
phylonode:
 - type_id ? (root,leaf,internal)
"""
@bcs_session(read_only=False)
def import_phylotrees(input_file, name = None, comment = None, analysis_id = None):
    # Check newick format
    from Bio import Phylo
    try:
        tree = Phylo.read(input_file, 'newick')
    except Exception as e:
        return f'The file could not be loaded correctly. Check that it is the newick format.\n{e}', 500
    # Create the not null dbxref for this phylotree
    if not name:
        name = os.path.basename(input_file)
    # Create the new phylotree
    phylotree = __new_phylotree(name, comment, analysis_id)
    bcs_phylotree = __phylotree2bcs(phylotree)
    # Get phylonodes insertion
    phylonodes = __tree2phylonodes(phylotree.phylotree_id, tree.root, None, [0])
    # try:
    #     chado_session.commit()
    # except Exception as e:
    #     chado_session.rollback()
    #     print(e)
    #     return f'The phylotree could not be inserted correctly.\n{e}', 500
    return 'The phylotree was imported correctly.', 200

def __new_phylotree(name, comment = None, analysis_id = None):
    dbxref = __new_phylotree_dbxref(name)
    from biobarcoding.db_models.chado import Phylotree, Analysis
    phylotree = Phylotree(dbxref_id=dbxref.dbxref_id, name=name)
    if comment:
        phylotree.comment = comment
    if analysis_id:
        phylotree.analysis_id = chado_session.query(Analysis).filter(Analysis.analysis_id==analysis_id).first().analysis_id
    chado_session.add(phylotree)
    chado_session.flush()
    return phylotree

def __new_phylotree_dbxref(name):
    from biobarcoding.db_models.chado import Db, Dbxref
    dbxref = Dbxref(
        db_id=chado_session.query(Db).filter(Db.name=='null').first().db_id,
        accession=f'phylotree:{name}',
        version=time.strftime("%Y %b %d %H:%M:%S"))
    chado_session.add(dbxref)
    chado_session.flush()
    return dbxref

def __phylotree2bcs(phylotree):
    from biobarcoding.services import get_or_create
    from biobarcoding.db_models.bioinformatics import PhylogeneticTree
    bcs_phylotree = get_or_create(db_session, PhylogeneticTree,
        chado_phylotree_id = phylotree.phylotree_id,
        chado_table = 'phylotree',
        name = phylotree.name)
    db_session.merge(bcs_phylotree)
    # db_session.commit()
    return bcs_phylotree

def __tree2phylonodes(phylotree_id, node, parent_id=None, index=[0]):
    from biobarcoding.db_models.chado import Phylonode, Feature
    phylonodes = []
    feature_id = chado_session.query(Feature.feature_id).filter(Feature.uniquename==node.name).first()
    phylonode = Phylonode(
        phylotree_id=phylotree_id,
        parent_phylonode_id=parent_id,
        feature_id=feature_id,
        label=node.name,
        distance=node.branch_length,
        left_idx=index[0],
        right_idx=index[0]+1)
    index[0]+=1
    chado_session.add(phylonode)
    chado_session.flush()
    if len(node.clades)>0:
        for clade in node.clades:
            phylonodes += __tree2phylonodes(phylotree_id, clade, phylonode.phylonode_id, index)
        phylonode.right_idx = index[0]
        chado_session.add(phylonode)
        chado_session.flush()
    index[0]+=1
    return phylonodes + [phylonode]


def __get_query(phylotree_id = None, analysis_id = None, name = None, comment = None, feature_id = None):
    from biobarcoding.db_models.chado import Phylotree, Dbxref
    query = chado_session.query(Phylotree)\
        .join(Dbxref)\
        .filter(Dbxref.accession!='taxonomy')
    if phylotree_id:
        query = query.filter(Phylotree.phylotree_id==phylotree_id)
    if analysis_id:
        query = query.filter(Phylotree.analysis_id==analysis_id)
    if name:
        query = query.filter(Phylotree.name==name)
    if comment:
        query = query.filter(Phylotree.comment==comment)
    if feature_id:
        from biobarcoding.db_models.chado import Phylonode
        query = query.join(Phylonode)\
            .filter(Phylonode.feature_id==feature_id)
    return query
