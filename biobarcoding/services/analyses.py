from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.db_models.chado import Analysis
from biobarcoding.rest import Issue, IType, filter_parse, paginator


def create(**kwargs):
    content = None
    try:
        chado_session.add(Analysis(**kwargs))
        issues, status = [Issue(IType.INFO, f'CREATE analyses: The analysis "{kwargs.get("program")} {kwargs.get("programversion")}" created successfully.\{res}')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE analyses: The analysis "{kwargs.get("program")} {kwargs.get("programversion")}" could not be created.')], 500
    return issues, content, status


def read(id=None, **kwargs):
    content = None
    try:
        content = __get_query(id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ analyses: The analyses were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ analyses: The analyses could not be read.')], 500
    return issues, content, status


def update(id, **kwargs):
    issues = [Issue(IType.WARNING, 'UPDATE alignments: dummy completed')]
    return issues, None, 200


def delete(id=None, **kwargs):
    content = None
    try:
        query = __get_query(id, **kwargs)
        # from biobarcoding.services.sequences import delete as delete_sequences
        # _ids = [msa.analysis_id for msa in query.all()]
        # delete_sequences(filter={'analysis_id':{'op':'in','analysis_id':_ids}})
        # __delete_from_bcs(_ids)
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
    if filter.get('feature_id'):
        from biobarcoding.db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.analysis_id)\
            .filter(filter_parse(AnalysisFeature, {'feature_id':filter.get('feature_id')})).all()
        clause.append(Analysis.analysis_id.in_(_ids))
    return clause


def __get_query_ordered(query, order):
    # query = query.order(order_parse(Analysis, kwargs.get('order'), __aux_own_order))
    return query