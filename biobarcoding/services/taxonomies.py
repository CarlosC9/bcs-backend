from ..db_models import DBSessionChado as chado_session
from ..db_models.chado import Phylotree
from ..rest import Issue, IType, filter_parse
from . import get_query


def create(**kwargs):
    return {'status': 'success', 'message': 'CREATE: taxonomies dummy completed'}, 200


def read(id=None, **kwargs):
    content, count = None, 0
    try:
        content, count = __get_query(id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ taxonomies: The taxonomies were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.INFO, 'READ taxonomies: The taxonomies could not be read.')], 400
    return issues, content, count, status


def update(id, **kwargs):
    issues = [Issue(IType.WARNING, 'UPDATE taxonomies: dummy completed')]
    content = None
    return issues, content, 200


def delete(id=None, **kwargs):
    content = None
    try:
        content, count = __get_query(id, **kwargs)
        content = content.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO,
                                f'DELETE taxonomies: The {content} taxonomies were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'DELETE taxonomies: The taxonomies could not be removed.')], 404
    return issues, content, status


def import_file(input_file, format='obo', **kwargs):
    content = None
    try:
        from flask import current_app
        cfg = current_app.config
        named = f' -n {kwargs.get("name")} ' if kwargs.get("name") else ''
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
        chado_session.execute(
            "SELECT setval('phylonode_phylonode_id_seq', (SELECT MAX(phylonode_id) FROM phylonode)+1);")
        # chado_session.execute("ALTER SEQUENCE phylonode_phylonode_id_seq RESTART WITH (SELECT MAX(phylonode_id) FROM phylonode)+1;")
        issues, status = [Issue(IType.INFO,
                                f'IMPORT taxonomies: The taxonomy {os.path.basename(input_file)} was successfully imported.')], 200
        if err:
            issues, status = [Issue(IType.INFO,
                                    f'IMPORT taxonomies: The taxonomy {os.path.basename(input_file)} was barely imported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR,
                                f'IMPORT taxonomies: The taxonomy {os.path.basename(input_file)} could not be imported.')], 409
    return issues, content, status


def export(id=None, **kwargs):
    issues = [Issue(IType.WARNING, 'EXPORT taxonomies: dummy completed')]
    content = None
    return issues, content, 200


def __get_query(phylotree_id=None, **kwargs):
    if phylotree_id:
        query = chado_session.query(Phylotree).filter(Phylotree.phylotree_id == phylotree_id)
        return query, query.count()
    from biobarcoding.db_models.chado import Dbxref
    dbxref_id = chado_session.query(Dbxref.dbxref_id).filter(Dbxref.accession == 'taxonomy').first()
    return get_query(chado_session, Phylotree, **kwargs, dbxref_id=dbxref_id,
                     aux_filter=__aux_own_filter, aux_order=__aux_own_order)


def __aux_own_filter(filter):
    clause = []
    # if 'analysis_id' in filter:
    #     from biobarcoding.db_models.chado import AnalysisFeature
    #     _ids = chado_session.query(AnalysisFeature.feature_id)\
    #             .filter(
    #                 filter_parse(AnalysisFeature, [{'analysis_id': filter.get('analysis_id')}]))
    #     clause.append(Phylotree.phylotree_id.in_(_ids))
    return clause


def __aux_own_order(order):
    # query = query.order(order_parse(Feature, kwargs.get('order'), __aux_own_order))
    return []
