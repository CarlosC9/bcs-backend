from ..db_models import DBSession as db_session
from ..db_models.sysadmin import BrowserFilter
from ..rest import Issue, IType
from . import get_query


def create(datatype, **kwargs):
    try:
        if not kwargs.get('name') or not datatype:
            raise Exception
        # get Identity
        content = BrowserFilter(**kwargs, type=datatype)
        from flask import g
        content.user_id = g.n_session.identity.id
        db_session.add(content)
        issues, status = [Issue(IType.INFO, f'CREATE browser_filters: It was created successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE browser_filters: It could not be created.')], 409
    return issues, content, status


def read(datatype, id=None, **kwargs):
    content, count = None, 0
    try:
        content, count = __get_query(datatype, id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ browser_filters: The browser_filters were read successfully.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ browser_filters: The browser_filters could not be read.')], 400
    return issues, content, count, status


def update(datatype, id, **kwargs):
    content = None
    try:
        content, count = __get_query(datatype, id)
        content = content.update(kwargs)
        issues, status = [Issue(IType.INFO, f'UPDATE browser_filters: It was successfully updated.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'UPDATE browser_filters: It could not be updated.')], 409
    return issues, content, status


def delete(datatype, id=None, **kwargs):
    content = None
    try:
        content, count = __get_query(datatype, id, **kwargs)
        content = content.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO,
                                f'DELETE browser_filters: The {content} browser_filters were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE browser_filters: The browser_filters could not be removed.')], 404
    return issues, content, status


def __get_query(type=None, id=None, **kwargs):
    return get_query(db_session, BrowserFilter, id=id, type=type, **kwargs,
                     aux_filter=__aux_own_filter, aux_order=__aux_own_order)


def __aux_own_filter(filter):
    return []


def __aux_own_order(order):
    return []


def read_form(datatype):
    content = None
    try:
        from biobarcoding.forms.filter_forms import getFilterSchema
        content = getFilterSchema(datatype)
        issues, status = [Issue(IType.INFO, 'READ browser_filter forms: successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ browser_filter forms: could not be read.')], 400
    return issues, content, status
