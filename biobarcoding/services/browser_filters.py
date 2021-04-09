from biobarcoding.db_models import DBSession as bcs_session
from biobarcoding.db_models.sysadmin import BrowserFilter
from biobarcoding.rest import Issue, IType, filter_parse, paginator


def create(datatype, **kwargs):
    try:
        if not kwargs.get('name') or not datatype:
            raise Exception
        # get Identity
        new_filter = BrowserFilter(**kwargs, type=datatype)
        from flask import g
        new_filter.user_id = g.bcs_session.identity.id
        bcs_session.add(new_filter)
        issues, status = [Issue(IType.INFO, f'CREATE browser_filters: It was created successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE browser_filters: It could not be created.')], 500
    return issues, None, status


count = 0
def read(datatype, id=None, **kwargs):
    content = None
    try:
        content = __get_query(datatype, id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ browser_filters: The browser_filters were read successfully.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ browser_filters: The browser_filters could not be read.')], 500
    return issues, content, count, status


def update(datatype, id, **kwargs):
    try:
        __get_query(datatype, id).update(kwargs)
        issues, status = [Issue(IType.INFO, f'UPDATE browser_filters: It was updated successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'UPDATE browser_filters: It could not be updated.')], 500
    return issues, None, status


def delete(datatype, id=None, **kwargs):
    content = None
    try:
        resp = __get_query(datatype, id, **kwargs).delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE browser_filters: The {resp} browser_filters were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE browser_filters: The browser_filters could not be removed.')], 500
    return issues, content, status


def __get_query(type=None, id=None, **kwargs):
    query = bcs_session.query(BrowserFilter)
    global count
    count = 0
    if id:
        query = query.filter(BrowserFilter.id == id)
    else:
        query = query.filter(BrowserFilter.type == type)
        if 'filter' in kwargs:
            query = query.filter(filter_parse(BrowserFilter, kwargs.get('filter'), __aux_own_filter))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            count = query.count()
            query = paginator(query, kwargs.get('pagination'))
    return query


def __aux_own_filter(filter):
    clause = []
    # if filter.get('feature_id'):
    #     from biobarcoding.db_models.chado import AnalysisFeature
    #     _ids = bcs_session.query(AnalysisFeature.analysis_id)\
    #         .filter(filter_parse(AnalysisFeature, {'feature_id':filter.get('feature_id')})).all()
    #     clause.append(BrowserFilter.analysis_id.in_(_ids))
    return clause


def __get_query_ordered(query, order):
    # query = query.order(order_parse(BrowserFilter, kwargs.get('order'), __aux_own_order))
    return query


def read_form(datatype):
    content = None
    try:
        from biobarcoding.forms.filter_forms import getFilterSchema
        content = getFilterSchema(datatype)
        issues, status = [Issue(IType.INFO, 'READ browser_filter forms: successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ browser_filter forms: could not be read.')], 500
    return issues, content, status