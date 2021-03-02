from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.db_models.chado import Analysis
from biobarcoding.rest import Issue, IType, filter_parse, paginator


def create_analyses(program, programversion, name = None, sourcename = None, description = None, algorithm = None, sourceversion = None, sourceuri = None, timeexecuted = None):
    content = { 'program':program, 'programversion':programversion, 'name':name, 'sourcename':sourcename, 'description':description, 'algorithm':algorithm,
                      'sourceversion':sourceversion, 'sourceuri':sourceuri, 'timeexecuted':timeexecuted }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        from biobarcoding.services import conn_chado
        conn = conn_chado()
        res = conn.analysis.add_analysis(
            program = program,
            programversion = programversion,
            name = name,
            description = description,
            algorithm = algorithm,
            sourcename = sourcename,
            sourceversion = sourceversion,
            sourceuri = sourceuri,
            date_executed = timeexecuted)
        issues, status = [Issue(IType.INFO, f'CREATE analyses: The analysis "{program} {programversion}" created successfully.\{res}')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE analyses: The analysis "{program} {programversion}" could not be created.')], 500
    return issues, content, status


def read_analyses(analysis_id=None, ids=None, name=None, program=None, programversion=None, algorithm=None,
                  sourcename=None, sourceversion=None, sourceuri=None, description=None, feature_id=None):
    content = { 'analysis_id':analysis_id, 'ids':ids, 'name':name, 'program':program, 'programversion':programversion,
                'algorithm':algorithm, 'sourcename':sourcename, 'sourceversion':sourceversion, 'sourceuri':sourceuri,
                'description':description, 'feature_id':feature_id}
    content = {k:v for k,v in content.items() if v is not None}
    try:
        content = __get_query(analysis_id=analysis_id, ids=ids, name=name, program=program, programversion=programversion,
            algorithm=algorithm, sourcename=sourcename, sourceversion=sourceversion, sourceuri=sourceuri,
            description=description, feature_id=feature_id)
        if analysis_id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ analyses: The analyses were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ analyses: The analyses could not be read.')], 500
    return issues, content, status


def update_analyses(analysis_id, program, programversion, name = None, description = None, algorithm = None, sourcename = None, sourceversion = None, sourceuri = None, timeexecuted = None):
    issues = [Issue(IType.WARNING, 'UPDATE alignments: dummy completed')]
    content = { 'analysis_id':analysis_id, 'program':program, 'programversion':programversion, 'name':name,
                'sourcename':sourcename, 'description':description, 'algorithm':algorithm,
                'sourceversion':sourceversion, 'sourceuri':sourceuri, 'timeexecuted':timeexecuted }
    return issues, content, 200


def delete_analyses(analysis_id=None, ids=None, name=None, program=None, programversion=None, algorithm=None, sourcename=None, sourceversion=None, sourceuri=None, description=None):
    content = { 'analysis_id':analysis_id, 'ids':' '.join(ids) if ids else None, 'name':name, 'program':program,
                'programversion':programversion, 'algorithm':algorithm, 'sourcename':sourcename,
                'sourceversion':sourceversion, 'sourceuri':sourceuri, 'description':description }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        query = __get_query(analysis_id=analysis_id, ids=ids, name=name, program=program, programversion=programversion,
                          algorithm=algorithm, sourcename=sourcename, sourceversion=sourceversion, sourceuri=sourceuri,
                          description=description)
        # from biobarcoding.services.sequences import delete_sequences
        # for msa in query.all():
        #     delete_sequences(analysis_id=msa.analysis_id)
        #     __delete_from_bcs(msa.analysis_id)
        resp = query.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE alignments: The {resp} alignments were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE alignments: The alignments could not be removed.')], 500
    return issues, content, status


def __get_query(id=None, **kwargs):
    query = chado_session.query(Analysis)
    if id:
        query = query.filter(Analysis.analysis_id == id)
    else:
        if 'filter' in kwargs:
            query = query.filter(filter_parse(Analysis, kwargs.get('filter'), __aux_own_filter))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            query = paginator(query, kwargs.get('pagination'))
    return query


def __aux_own_filter(filter):
    clause = []
    if 'feature_id' in filter and filter.get('feature_id'):
        from biobarcoding.db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.analysis_id)\
            .filter(AnalysisFeature.feature_id==filter.get('feature_id'))
        clause.append(Analysis.analysis_id.in_(_ids))
    return clause


def __get_query_ordered(query, order):
    # query = query.order(order_parse(Analysis, kwargs.get('order'), __aux_own_order))
    return query