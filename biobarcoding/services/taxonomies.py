from biobarcoding.db_models import DBSessionChado as chado_session

def create_taxonomies(name, comment = None):
    return {'status':'success','message':'CREATE: taxonomies dummy completed'}, 200


def read_taxonomies(id = None):
    from biobarcoding.services import chado2json
    from biobarcoding.db_models.chado import Phylotree, Dbxref
    result = chado_session.query(Phylotree)\
        .join(Dbxref)\
        .filter(Dbxref.accession=='taxonomy')
    if id:
        result = result.filter(Phylotree.phylotree_id==id)
        return chado2json(result)[0], 200
    return chado2json(result), 200


def update_taxonomies(id, name = None, comment = None):
    return {'status':'success','message':'UPDATE: taxonomies dummy completed'}, 200


def delete_taxonomies(id=None, ids=None):
    from biobarcoding.db_models.chado import Phylotree, Dbxref
    result = chado_session.query(Phylotree)\
        .filter(Phylotree.dbxref_id==chado_session.query(Dbxref.dbxref_id)\
                .filter(Dbxref.accession=='taxonomy').first())
    msg='[ '
    if id:
        result.filter(Phylotree.phylotree_id==id)
        msg+=f'{id} '
    if ids:
        result.filter(Phylotree.phylotree_id.in_(ids))
        msg+=f'{ids} '
    msg+=']'
    result.delete()
    chado_session.commit()
    return {'status':'success','message':f'Taxonomies deleted: {msg}'}, 200


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
    chado_session.execute("SELECT setval('phylonode_phylonode_id_seq', (SELECT MAX(phylonode_id) FROM phylonode)+1);")
    # chado_session.execute("ALTER SEQUENCE phylonode_phylonode_id_seq RESTART WITH (SELECT MAX(phylonode_id) FROM phylonode)+1;")
    # chado_session.commit()
    if err:
        import os
        return {'status':'failure','message':f'Taxonomy in {os.path.basename(input_file)} could not be imported.\n{err}'}, 500
    return {'status':'success','message':f'Taxonomy in {os.path.basename(input_file)} imported properly.\n{out}'}, 200


def export_taxonomies(id = None, ids = None):
    return {'status':'success','message':'EXPORT: taxonomies dummy completed'}, 200
