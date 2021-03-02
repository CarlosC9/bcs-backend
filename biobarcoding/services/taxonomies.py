from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.rest import Issue, IType


def create_taxonomies(name, comment = None):
    return {'status':'success','message':'CREATE: taxonomies dummy completed'}, 200


def read_taxonomies(id = None, ids = None, name = None, comment = None):
    content = { 'id':id, 'ids':ids, 'name':name, 'comment':comment }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        query = __get_query()
        if id:
            content = query.first()
        else:
            content = query.all()
        issues, status = [Issue(IType.INFO, 'READ taxonomies: The taxonomies were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.INFO, 'READ taxonomies: The taxonomies could not be read.')], 500
    return issues, content, status


def update_taxonomies(id, name = None, comment = None):
    issues = [Issue(IType.WARNING, 'UPDATE taxonomies: dummy completed')]
    content = { 'id':id, 'name':name, 'comment':comment}
    return issues, content, 200


def delete_taxonomies(id=None, ids=None):
    content = { 'id':id, 'ids':ids }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        resp = __get_query(id=id, ids=ids).delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE taxonomies: The {resp} taxonomies were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'DELETE taxonomies: The taxonomies could not be removed.')], 500
    return issues, content, status


def import_taxonomies(input_file, name = None, comment = None):
    content = None
    try:
        from flask import current_app
        cfg = current_app.config
        named=''
        if name:
            named = f' -n {name} '
        import os
        dir_path = os.path.dirname(os.path.realpath(__file__))
        from biobarcoding.services import exec_cmds
        content, err = exec_cmds(f'''(cd {dir_path}/perl_scripts/ &&
            perl ./load_ncbi_taxonomy.pl\
                -H {cfg["CHADO_HOST"]}\
                -D {cfg["CHADO_DATABASE"]}\
                -u {cfg["CHADO_USER"]}\
                -p {cfg["CHADO_PASSWORD"]}\
                -d Pg\
                -i {input_file}\
                {named})''')
        chado_session.execute("SELECT setval('phylonode_phylonode_id_seq', (SELECT MAX(phylonode_id) FROM phylonode)+1);")
        # chado_session.execute("ALTER SEQUENCE phylonode_phylonode_id_seq RESTART WITH (SELECT MAX(phylonode_id) FROM phylonode)+1;")
        issues, status = Issue(IType.INFO, f'IMPORT taxonomies: The taxonomy {os.path.basename(input_file)} was successfully imported.'), 200
        if err:
            raise err
    except Exception as e:
        print(e)
        issues, status = Issue(IType.ERROR, f'IMPORT taxonomies: The taxonomy {os.path.basename(input_file)} could not be imported.'), 500
    return issues, content, status


def export_taxonomies(id = None, ids = None, format = None):
    issues = [Issue(IType.WARNING, 'EXPORT taxonomies: dummy completed')]
    content = { 'id':id, 'ids':ids }
    return issues, content, 200


def __get_query(id = None, ids = None, name = None, comment = None):
    from biobarcoding.db_models.chado import Phylotree, Dbxref
    result = chado_session.query(Phylotree)\
        .join(Dbxref)\
        .filter(Dbxref.accession=='taxonomy')
    if id:
        result = result.filter(Phylotree.phylotree_id==id)
    if ids:
        result = result.filter(Phylotree.phylotree_id.in_(ids))
    if name:
        result = result.filter(Phylotree.name==name)
    if comment:
        result = result.filter(Phylotree.comment==comment)
    return result
