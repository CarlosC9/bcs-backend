import os, time
from biobarcoding.db_models import DBSession as bcs_session
from biobarcoding.db_models import DBSessionChado as chado_session

def create_phylotrees(name = None, comment = None):
    return {'status':'success','message':'CREATE: phylotrees dummy completed'}, 200

def read_phylotrees(id = None):
    from biobarcoding.db_models.chado import Phylotree, Dbxref
    result = chado_session.query(Phylotree)\
        .join(Dbxref)\
        .filter(Dbxref.accession!='taxonomy')
    if id:
        result = result.filter(Phylotree.phylotree_id==id)
    response = []
    for value in result.all():
        tmp = value.__dict__
        tmp.pop('_sa_instance_state', None)
        response.append(tmp)
    if id:
        return response[0], 200
    return response, 200

def update_phylotrees(phylotree_id, name = None, comment = None, analysis_id = None):
    return {'status':'success','message':'UPDATE: phylotrees dummy completed'}, 200


def delete_phylotrees(phylotree_id = None):
    return {'status':'success','message':'DELETE: phylotrees dummy completed'}, 200

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


def export_phylotrees(phylotree_id = None):
    return {'status':'success','message':'UPDATE: phylotrees dummy completed'}, 200

"""
# NGD newick phylotree import
__TODO__
phylotree:
 - type_id ? cvterm['phylogeny'].cvterms_id
 - for each selected cvterms: new phylotreeprop
 - store the input_file content into a phylotreeprop EDAM_newick
phylonode:
 - type_id ? (root,leaf,internal)
"""
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
    try:
        chado_session.commit()
    except Exception as e:
        chado_session.rollback()
        print(e)
        return f'The phylotree could not be inserted correctly.\n{e}', 500
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
    bcs_phylotree = get_or_create(bcs_session, PhylogeneticTree,
        chado_phylotree_id = phylotree.phylotree_id,
        chado_table = 'phylotree',
        bcs_phylotree.name = phylotree.name)
    bcs_session.merge(bcs_phylotree)
    bcs_session.commit()
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
