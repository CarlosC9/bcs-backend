def create_phylotree(input_file, name = None, comment = None, analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    # sequence_type='polypeptide'
    if not analysis_id:
        analysis_id = conn.analysis.get_analyses(name='Unknown analysis')[0]['analysis_id']
    res = conn.phylogeny.load_tree(input_file, analysis_id)
    return res

# sqlalchemy: phylotree
def read_phylotree(phylotree_id = None):
    pass

def update_phylotree(phylotree_id, name = None, comment = None, analysis_id = None):
    pass

def delete_phylotree(phylotree_id = None):
    pass
