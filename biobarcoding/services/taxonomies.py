def create_taxonomies(name, comment = None):
    return {'status':'success','message':'CREATE: taxonomies dummy completed'}, 200


def read_taxonomies(id = None):
    from biobarcoding.db_models import DBSessionChado
    from biobarcoding.db_models.chado import Phylotree, Dbxref
    result = DBSessionChado().query(Phylotree)\
        .join(Dbxref)\
        .filter(Dbxref.accession=='taxonomy')
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


def update_taxonomies(id, name = None, comment = None):
    return {'status':'success','message':'UPDATE: taxonomies dummy completed'}, 200


def delete_taxonomies(id):
    return {'status':'success','message':'DELETE: taxonomies dummy completed'}, 200


def import_taxonomies(input_file, name = None, comment = None):
    from flask import current_app
    cfg = current_app.config
    named=''
    if name:
        named = f' -n {name} '
    from biobarcoding.services import exec_cmds
    out, err = exec_cmds([f'''(cd ./biobarcoding/services/perl_scripts/ &&
        perl ./load_ncbi_taxonomy.pl\
            -H {cfg["CHADO_HOST"]}\
            -D {cfg["CHADO_DATABASE"]}\
            -u {cfg["CHADO_USER"]}\
            -p {cfg["CHADO_PASSWORD"]}\
            -d Pg\
            -i {input_file}\
            {named})'''])
    if err:
        import os
        return {'status':'failure','message':f'Taxonomy in {os.path.basename(input_file)} could not be imported.\n{err}'}, 500
    return {'status':'success','message':f'Taxonomy in {os.path.basename(input_file)} imported properly.\n{out}'}, 200


def export_taxonomies(id = None):
    return {'status':'success','message':'EXPORT: taxonomies dummy completed'}, 200
