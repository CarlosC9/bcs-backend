from ... import get_or_create, get_query, get_orm_params
from ....db_models import DBSessionChado as chado_session
from ....db_models.chado import Analysis, AnalysisCvterm
from ....rest import Issue, IType, filter_parse


##
# CREATE
##

def __check_ansis_values(**values):
    if values.get('job_id'):
        values['sourcename'], values['sourceversion'], values['sourceuri'] = values.get('job_id'), 'job', f'/jobs/{values.get("job_id")}'
    if not (values.get('program') or values.get('programversion') or values.get('sourcename')):
        raise Exception('Missing required params ("program", "programversion", "sourcename").')
    if not values.get('program'):
        values['program'] = 'unknown'
    if not values.get('programversion'):
        values['programversion'] = 'unknown'
    if not values.get('name'):
        values['name'] = f"{values['program']} {values['programversion']}"
    if not values.get('sourcename'):
        values['sourcename'] = 'unknown'
    return get_orm_params(Analysis, **values)


def create(**kwargs):
    content = None
    try:
        values = __check_ansis_values(**kwargs)
        content = Analysis(**values)
        chado_session.add(content)
        issues, status = [Issue(IType.INFO,
                                f'CREATE analyses: The analysis "{kwargs.get("program")} {kwargs.get("programversion")}" created successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR,
                                f'CREATE analyses: The analysis "{kwargs.get("program")} {kwargs.get("programversion")}" could not be created.')], 409
    return issues, content, status


##
# READ
##

def read(id=None, **kwargs):
    content, count = None, 0
    try:
        content, count = __get_query(id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ analyses: The analyses were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ analyses: The analyses could not be read.')], 400
    return issues, content, count, status


##
# UPDATE
##

def update(id, **kwargs):
    content = None
    try:
        content = __get_query(id)[0].one().update(kwargs)
        issues, status = [Issue(IType.INFO, f'UPDATE analyses: The analysis "{id}" was successfully updated.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'UPDATE analyses: The analysis "{id}" could not be updated.')], 409
    return issues, content, status


##
# DELETE
##

def delete(id=None, **kwargs):
    content = None
    try:
        query, count = __get_query(id, **kwargs)
        # from ..bos.sequences import delete as delete_sequences
        # _ids = [msa.analysis_id for msa in query.all()]
        # delete_sequences(filter={'analysis_id':{'op':'in','analysis_id':_ids}})
        # __delete_from_bcs(_ids)
        resp = query.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE analyses: The {resp} analyses were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE analyses: The analyses could not be removed.')], 404
    return issues, content, status


##
# GETTER AND OTHERS
##

def __get_query(analysis_id=None, **kwargs):
    if analysis_id:
        query = chado_session.query(Analysis).filter(Analysis.analysis_id == analysis_id)
        return query, query.count()
    if kwargs.get('job_id'):
        # TODO: will there be multiple analyzes for a single job ?
        query = chado_session.query(Analysis).filter(Analysis.sourcename == str(kwargs.get('job_id')),
                                                     Analysis.sourceversion == 'job')
        return query, query.count()
    return get_query(chado_session, Analysis, aux_filter=__aux_ansis_filter, aux_order=__aux_ansis_order, **kwargs)


def __aux_ansis_filter(filter):
    clause = []

    if filter.get('job_id'):
        from ....db_models.chado import Analysis
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, {'sourcename': filter.get('job_id'),
                                            'sourceversion': {'op': 'eq', 'unary': 'job'}}))
        clause.append(Analysis.analysis_id.in_(_ids))

    if filter.get('feature_id'):
        from ....db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.analysis_id) \
            .filter(filter_parse(AnalysisFeature, {'feature_id': filter.get('feature_id')}))
        clause.append(Analysis.analysis_id.in_(_ids))

    if filter.get('organism_id'):
        from ....db_models.chado import Feature
        _ids = chado_session.query(Feature.feature_id) \
            .filter(filter_parse(Feature, {'organism_id': filter.get('organism_id')}))
        from ....db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.analysis_id) \
            .filter(AnalysisFeature.feature_id.in_(_ids))
        clause.append(Analysis.analysis_id.in_(_ids))

    if 'phylotree_id' in filter:
        from ....db_models.chado import Phylotree
        _ids = chado_session.query(Phylotree.analysis_id) \
            .filter(filter_parse(Phylotree, [{'phylotree_id': filter.get('phylotree_id')}]))
        clause.append(Analysis.analysis_id.in_(_ids))

    if "cvterm_id" in filter:
        from ....db_models.chado import AnalysisCvterm
        _ids = chado_session.query(AnalysisCvterm.analysis_id) \
            .filter(filter_parse(AnalysisCvterm, [{'cvterm_id': filter.get('cvterm_id')}]))
        clause.append(Analysis.analysis_id.in_(_ids))

    if "prop_cvterm_id" in filter:
        from ....db_models.chado import Analysisprop
        _ids = chado_session.query(Analysisprop.analysis_id) \
            .filter(filter_parse(Analysisprop, [{'type_id': filter.get('prop_cvterm_id')}]))
        clause.append(Analysis.analysis_id.in_(_ids))

    from datetime import datetime
    if "added-from" in filter:
        filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, {'timeexecuted': filter.get("added-from")}))
        clause.append(Analysis.analysis_id.in_(_ids))
    if "added-to" in filter:
        filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, {'timeexecuted': filter.get("added-to")}))
        clause.append(Analysis.analysis_id.in_(_ids))

    return clause


def __aux_ansis_order(order):
    clauses = []
    return clauses
