def create_phylotrees(input_file, name = None, comment = None, analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    # sequence_type='polypeptide'
    if not analysis_id:
        analysis_id = conn.analysis.get_analyses(name='Unknown analysis')[0]['analysis_id']
    res = conn.phylogeny.load_tree(input_file, analysis_id)
    return res

def read_phylotrees(id = None):
    from biobarcoding.db_models import DBSessionChado
    from biobarcoding.db_models.chado import Phylotree, Cvterm
    result = DBSessionChado().query(Phylotree)
    if id:
        result = result.filter(Phylotree.phylotree_id==id)
    response = []
    for value in result.all():
        tmp = value.__dict__
        tmp.pop('_sa_instance_state', None)
        response.append(tmp)
    return {'status':'success','message':response}, 200

def update_phylotrees(phylotree_id, name = None, comment = None, analysis_id = None):
    pass

def delete_phylotrees(phylotree_id = None):
    pass
