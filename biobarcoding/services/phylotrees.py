import os, time

from biobarcoding.db_models import DBSession as db_session
from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.rest import Issue, IType
from biobarcoding.db_models.chado import Phylotree, Phylonode

def create_phylotrees(name = None, comment = None):
    issues = [Issue(IType.WARNING, 'CREATE phylotrees: dummy completed')]
    content = { 'name':name, 'comment':comment }
    return issues, content, 200


def read_phylotrees(id = None, analysis_id = None, name = None, comment = None, feature_id = None):
    content = { 'id':id, 'analysis_id':analysis_id, 'name':name, 'comment':comment, 'feature_id':feature_id }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        content = __get_query(id, analysis_id=analysis_id, name=name, comment=comment, feature_id=feature_id)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ phylotrees: The phylotrees were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ phylotrees: The phylotrees could not be read.')], 500
    return issues, content, status



def update_phylotrees(phylotree_id, name = None, comment = None, analysis_id = None):
    issues = [Issue(IType.WARNING, 'UPDATE phylotrees: dummy completed')]
    content = { 'phylotree_id':phylotree_id, 'name':name, 'comment':comment, 'analysis_id':analysis_id }
    return issues, content, 200


def __delete_from_bcs(*args):
    from biobarcoding.db_models.bioinformatics import PhylogeneticTree
    db_session.query(PhylogeneticTree).filter(PhylogeneticTree.chado_phylotree_id.in_(args))\
        .delete(synchronize_session='fetch')


def delete_phylotrees(id=None, ids=None, name = None, comment = None, analysis_id = None):
    content = { 'id':id, 'ids':ids, 'name':name, 'comment':comment, 'analysis_id':analysis_id}
    content = {k:v for k,v in content.items() if v is not None}
    try:
        query = __get_query(phylotree_id=id, ids=ids)
        for phylo in query.all():
            __delete_from_bcs(phylo.phylotree_id)
        resp = query.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE phylotrees: The {resp} phylotrees were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE phylotrees: The phylotrees could not be removed.')], 500
    return issues, content, status


# REMOVE Legacy import with python-chado
# def __import_phylotrees(input_file, name = None, comment = None, analysis_id = None):
#     from biobarcoding.services import conn_chado
#     conn = conn_chado()
#     try:
#         if not analysis_id:
#             analysis_id = conn.analysis.get_analyses()[0]['analysis_id']
#         response = conn.phylogeny.load_tree(input_file, analysis_id, name=name)
#         return {'status':'success','message':f'{response} phylotrees were successfully imported.'}, 200
#     except Exception as e:
#         print(e)
#         return {'status':'failure','message':'The phylotree could not be imported.'}, 500


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
def import_phylotrees(input_file, name=None, comment=None, analysis_id=None, format='newick'):
    content = { 'input_file':input_file, 'name':name, 'comment':comment, 'analysis_id':analysis_id }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        # Check newick format
        from Bio import Phylo
        tree = Phylo.read(input_file, format)
    except Exception as e:
        issues = [Issue(IType.ERROR, f'IMPORT phylotress: The file {input_file}.{format} could not be imported.')]
        return issues, content, 500
    try:
        if not name:
            name = os.path.basename(input_file)
        # Create the new phylotree
        phylotree = __new_phylotree(name, comment, analysis_id)
        bcs_phylotree = __phylotree2bcs(phylotree)
        # Get phylonodes insertion
        phylonodes = __tree2phylonodes(phylotree.phylotree_id, tree.root, None, [0])
        issues, status = [Issue(IType.INFO, 'IMPORT phylotrees: The phylotree was successfully imported.')], 200
    except Exception as e:
        issues, status = [Issue(IType.ERROR, f'IMPORT phylotress: The file {input_file}.{format} could not be imported.')], 500
    return issues, content, status

def __new_phylotree(name, comment = None, analysis_id = None):
    dbxref = __new_phylotree_dbxref(name)
    from biobarcoding.db_models.chado import Analysis
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


# NGD newick phylotree export
# TODO ¿rebuilt phylogeny from chado?
def __tree2newick(node, is_root=True):
    result=''
    children = chado_session.query(Phylonode).filter(Phylonode.parent_phylonode_id==node.phylonode_id)
    if children.count()>0:
        result+='(' if is_root else '\n('
        i=1
        for n in children.all():
            result+=__tree2newick(n, False)
            result+=',' if i<children.count() else ')'
            i+=1
    if node.label or node.distance:
        label = node.label if node.label else ''
        distance = node.distance if node.distance else '0.00000'
        result+=f'\n{label}:{distance}'
    result+=';' if is_root else ''
    return result


def __tree2file(phylotree_id, format, output_file):
    result = ''
    if format=='newick':
        root = chado_session.query(Phylonode).filter(Phylonode.phylotree_id==phylotree_id,
                                                     Phylonode.parent_phylonode_id==None).first()
        result=__tree2newick(root)
    with open(output_file, "w") as file:
        file.write(result)
    return output_file


def export_phylotrees(phylotree_id=None, format='newick', output_file='/tmp/output_seqs.fas'):
    try:
        if __get_query(phylotree_id).first():
            __tree2file(phylotree_id, format, output_file)
        issues, status = [Issue(IType.INFO, 'EXPORT phylotrees: The phylotree were successfully exported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'EXPORT phylotrees: The phylotree could not be exported.')], 500
    return issues, output_file, status


def __get_query(phylotree_id = None, ids=None, analysis_id = None, name = None, comment = None, feature_id = None):
    from biobarcoding.db_models.chado import Dbxref
    dbxref_id = chado_session.query(Dbxref.dbxref_id).filter(Dbxref.accession=='taxonomy').first()
    query = chado_session.query(Phylotree).filter(Phylotree.dbxref_id!=dbxref_id)
    if phylotree_id:
        query = query.filter(Phylotree.phylotree_id==phylotree_id)
    if ids:
        query = query.filter(Phylotree.phylotree_id.in_(ids))
    if analysis_id:
        query = query.filter(Phylotree.analysis_id==analysis_id)
    if name:
        query = query.filter(Phylotree.name==name)
    if comment:
        query = query.filter(Phylotree.comment==comment)
    if feature_id:
        query = query.join(Phylonode)\
            .filter(Phylonode.feature_id==feature_id)
    return query
