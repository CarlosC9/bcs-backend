from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.db_models.chado import Analysis
from biobarcoding.rest import Issue, IType, filter_parse, paginator


def create(**kwargs):
    content = None
    try:
        if not kwargs.get('name'):
            kwargs['name']= f"{kwargs['program']} {kwargs['programversion']}"
        chado_session.add(Analysis(**kwargs))
        issues, status = [Issue(IType.INFO, f'CREATE analyses: The analysis "{kwargs.get("program")} {kwargs.get("programversion")}" created successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE analyses: The analysis "{kwargs.get("program")} {kwargs.get("programversion")}" could not be created.')], 500
    return issues, content, status


count = 0
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
    return issues, content, count, status


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
    global count
    count = 0
    if id:
        query = query.filter(Analysis.analysis_id == id)
    else:
        if 'filter' in kwargs:
            query = query.filter(filter_parse(Analysis, kwargs.get('filter'), __aux_own_filter))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            count = query.count()
            query = paginator(query, kwargs.get('pagination'))
    return query


def __aux_own_filter(filter):
    clause = []

    if filter.get('feature_id'):
        from biobarcoding.db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.analysis_id)\
            .filter(filter_parse(AnalysisFeature, {'feature_id': filter.get('feature_id')}))
        clause.append(Analysis.analysis_id.in_(_ids))

    if filter.get('organism_id'):
        from biobarcoding.db_models.chado import Feature
        _ids = chado_session.query(Feature.feature_id)\
            .filter(filter_parse(Feature, {'organism_id': filter.get('organism_id')}))
        from biobarcoding.db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.analysis_id)\
            .filter(AnalysisFeature.feature_id.in_(_ids))
        clause.append(Analysis.analysis_id.in_(_ids))

    if 'phylotree_id' in filter:
        from biobarcoding.db_models.chado import Phylotree
        _ids = chado_session.query(Phylotree.analysis_id) \
            .filter(filter_parse(Phylotree, [{'phylotree_id': filter.get('phylotree_id')}]))
        clause.append(Analysis.analysis_id.in_(_ids))

    if "cvterm_id" in filter:
        from biobarcoding.db_models.chado import AnalysisCvterm
        _ids = chado_session.query(AnalysisCvterm.analysis_id) \
            .filter(filter_parse(AnalysisCvterm, [{'cvterm_id': filter.get('cvterm_id')}]))
        clause.append(Analysis.analysis_id.in_(_ids))

    if "prop_cvterm_id" in filter:
        from biobarcoding.db_models.chado import Analysisprop
        _ids = chado_session.query(Analysisprop.analysis_id) \
            .filter(filter_parse(Analysisprop, [{'type_id': filter.get('prop_cvterm_id')}]))
        clause.append(Analysis.analysis_id.in_(_ids))

    from datetime import datetime
    if "added-from" in filter:
        filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, {'timeexecuted':filter.get("added-from")}))
        clause.append(Analysis.analysis_id.in_(_ids))
    if "added-to" in filter:
        filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, {'timeexecuted':filter.get("added-to")}))
        clause.append(Analysis.analysis_id.in_(_ids))
    return clause


def __get_query_ordered(query, order):
    # query = query.order(order_parse(Analysis, kwargs.get('order'), __aux_own_order))
    return query