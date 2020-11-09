def create_phylotrees(name = None, comment = None):
    return {'status':'success','message':'CREATE: phylotrees dummy completed'}, 200


def read_phylotrees(id = None):
    from biobarcoding.db_models import DBSessionChado
    from biobarcoding.db_models.chado import Phylotree, Dbxref
    result = DBSessionChado().query(Phylotree)\
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


def import_phylotrees(input_file, name = None, comment = None, analysis_id = None):
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
