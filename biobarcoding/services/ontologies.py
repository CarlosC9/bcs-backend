from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.db_models.chado import Cvterm
from biobarcoding.rest import Issue, IType, filter_parse, paginator


def create(**kwargs):
    issues = [Issue(IType.WARNING, 'CREATE ontologies: dummy completed')]
    return issues, None, 200


count = 0
def read(id=None, **kwargs):
    content = None
    try:
        content = __get_query(id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ ontologies: The ontologies were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ ontologies: The ontologies could not be read.')], 400
    return issues, content, count, status


def update(id, remote_url = None, input_file = None, **kwargs):
    issues = [Issue(IType.WARNING, 'UPDATE ontologies: dummy completed')]
    return issues, None, 200


def delete(id = None, **kwargs):
    issues = [Issue(IType.WARNING, 'DELETE ontologies: dummy completed')]
    return issues, None, 200


def import_file(input_file, format='obo', **kwargs):
    content = None
    try:
        from flask import current_app
        cfg = current_app.config
        # f"""go2fmt.pl -p obo_text -w xml {input_file} | \
        #     go-apply-xslt oboxml_to_chadoxml - > {input_file}.xml"""
        # cmd = f"""go2chadoxml {input_file} > /tmp/{input_file}.chado.xml;
        #     stag-storenode.pl -d 'dbi:Pg:dbname={cfg['database']};host={cfg['host']};port=cfg['port']'
        #     --user {cfg['user']} --password {cfg['password']} /tmp/{input_file}.chado.xml"""
        import pronto
        import os
        namespace = pronto.Ontology(input_file).metadata.default_namespace
        if namespace:
            onto_name = f'-c {namespace}'
        else:
            onto_name = f'-c {os.path.basename(input_file)}'
        from biobarcoding.services import exec_cmds
        out, err = exec_cmds([
            f'''perl ./biobarcoding/services/perl_scripts/gmod_load_cvterms.pl\
                -H {cfg["CHADO_HOST"]}\
                -D {cfg["CHADO_DATABASE"]}\
                -r {cfg["CHADO_USER"]}\
                -p {cfg["CHADO_PASSWORD"]}\
                -d Pg -s null -u\
                {input_file}''',
            f'''perl ./biobarcoding/services/perl_scripts/gmod_make_cvtermpath.pl\
                -H {cfg["CHADO_HOST"]}\
                -D {cfg["CHADO_DATABASE"]}\
                -u {cfg["CHADO_USER"]}\
                -p {cfg["CHADO_PASSWORD"]}\
                -d Pg {onto_name}'''])
        issues, status = [Issue(IType.INFO, f'IMPORT ontologies: The {format} ontology were successfully imported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'IMPORT ontologies: The file {input_file} could not be imported.')], 409
    return issues, content, status


def export(id=None, format='obo', **kwargs):
    issues = [Issue(IType.WARNING, 'EXPORT ontologies: dummy completed')]
    return issues, None, 200


def __get_query(id=None, **kwargs):
    from biobarcoding.db_models.chado import Cv
    query = chado_session.query(Cv)
    global count
    count = 0
    if id:
        query = query.filter(Cv.cv_id == id)
    else:
        if 'filter' in kwargs:
            query = query.filter(filter_parse(Cv, kwargs.get('filter'), __aux_own_filter))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            count = query.count()
            query = paginator(query, kwargs.get('pagination'))
    return query

def __aux_own_filter(**kwargs):
    return []

def __get_query_ordered(query, **kwargs):
    return query


# TODO:
#  featureprop, analysisprop, phylotreeprop ?
#  phylotree ?
def read_cvterms(cv_id=None, cvterm_id=None, **kwargs):
    content = None
    try:
        result = __get_cvterm(cv_id, cvterm_id, **kwargs)
        if cvterm_id:
            content = result.first()
        else:
            content = result.all()
        issues, status = [Issue(IType.ERROR, f'READ ontology terms: The cvterms were successfully read.')], 200
    except Exception as e:
        issues, status = [Issue(IType.ERROR, f'READ ontology terms: The cvterms could not be read.')], 400
    return issues, content, status


def __get_cvterm(cv_id=None, cvterm_id=None, **kwargs):
    query = chado_session.query(Cvterm)
    if cv_id:
        query = query.filter(Cvterm.cv_id==cv_id)
    if cvterm_id:
        query = query.filter(Cvterm.cvterm_id==cvterm_id)
    else:
        if 'filter' in kwargs:
            query = query.filter(filter_parse(Cvterm, kwargs.get('filter'), __aux_cvterms_filter))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            query = paginator(query, kwargs.get('pagination'))
    return query


def __aux_cvterms_filter(filter):
    clause = []
    if filter.get('feature_id'):
        from biobarcoding.db_models.chado import Feature, Featureprop, FeatureCvterm
        _ids = chado_session.query(Feature.type_id)\
            .filter(Feature.feature_id==filter.get('feature_id')).all()
        _ids += chado_session.query(Featureprop.type_id)\
            .filter(Featureprop.feature_id==filter.get('feature_id')).all()
        _ids += chado_session.query(FeatureCvterm.cvterm_id)\
            .filter(FeatureCvterm.feature_id==filter.get('feature_id')).all()
        clause.append(Cvterm.cvterm_id.in_(_ids))
    if filter.get('analysis_id'):
        from biobarcoding.db_models.chado import Analysisprop, AnalysisCvterm
        _ids = chado_session.query(Analysisprop.type_id)\
            .filter(Analysisprop.analysis_id==filter.get('analysis_id')).all()
        _ids += chado_session.query(AnalysisCvterm.cvterm_id)\
            .filter(AnalysisCvterm.analysis_id==filter.get('analysis_id')).all()
        clause.append(Cvterm.cvterm_id.in_(_ids))
    if filter.get('phylotree_id'):
        from biobarcoding.db_models.chado import Phylotree, Phylotreeprop
        _ids = chado_session.query(Phylotree.type_id)\
            .filter(Phylotree.phylotree_id==filter.get('phylotree_id')).all()
        _ids += chado_session.query(Phylotreeprop.type_id)\
            .filter(Phylotreeprop.phylotree_id==filter.get('phylotree_id')).all()
        clause.append(Cvterm.cvterm_id.in_(_ids))
    return clause


def __get_query_ordered(query, order):
    # query = query.order(order_parse(Cvterm, kwargs.get('order'), __aux_own_order))
    return query
