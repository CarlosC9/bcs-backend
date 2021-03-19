from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.rest import Issue, IType, filter_parse, paginator

from biobarcoding.db_models.chado import Phylotree


def create(**kwargs):
    return {'status':'success','message':'CREATE: taxonomies dummy completed'}, 200


count = 0
def read(id = None, **kwargs):
    content = None
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
    return issues, content, count, status


def update(id, **kwargs):
    issues = [Issue(IType.WARNING, 'UPDATE taxonomies: dummy completed')]
    content = None
    return issues, content, 200


def delete(id=None, **kwargs):
    content = None
    try:
        resp = __get_query(id=id, **kwargs).delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE taxonomies: The {resp} taxonomies were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'DELETE taxonomies: The taxonomies could not be removed.')], 500
    return issues, content, status


def import_file(input_file, format = 'obo', **kwargs):
    content = None
    try:
        from flask import current_app
        cfg = current_app.config
        named=''
        if kwargs.get("name"):
            named = f' -n {kwargs.get("name")} '
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
        issues, status = [Issue(IType.INFO, f'IMPORT taxonomies: The taxonomy {os.path.basename(input_file)} was successfully imported.')], 200
        if err:
            issues, status = [Issue(IType.INFO, f'IMPORT taxonomies: The taxonomy {os.path.basename(input_file)} was barely imported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'IMPORT taxonomies: The taxonomy {os.path.basename(input_file)} could not be imported.')], 500
    return issues, content, status


def export(id = None, **kwargs):
    issues = [Issue(IType.WARNING, 'EXPORT taxonomies: dummy completed')]
    content = None
    return issues, content, 200


def __get_query(id = None, **kwargs):
    from biobarcoding.db_models.chado import Dbxref
    _tax_tag = chado_session.query(Dbxref.dbxref_id)\
        .filter(Dbxref.accession=='taxonomy').first()
    query = chado_session.query(Phylotree)\
        .filter(Phylotree.dbxref_id==_tax_tag)
    global count
    count = 0
    if id:
        query = query.filter(Phylotree.phylotree_id==id)
    else:
        if 'filter' in kwargs:
            query = query.filter(filter_parse(Phylotree, kwargs.get('filter')))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            count = query.count()
            query = paginator(query, kwargs.get('pagination'))
    return query


def __aux_own_filter(filter):
    clause=[]
    # if 'analysis_id' in filter:
    #     from biobarcoding.db_models.chado import AnalysisFeature
    #     _ids = chado_session.query(AnalysisFeature.feature_id)\
    #             .filter(
    #                 filter_parse(AnalysisFeature, [{'analysis_id': filter.get('analysis_id')}]))
    #     clause.append(Phylotree.phylotree_id.in_(_ids))
    return clause


def __get_query_ordered(query, order):
    # query = query.order(order_parse(Feature, kwargs.get('order'), __aux_own_order))
    return query