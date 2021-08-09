from ..db_models import DBSessionChado as chado_session
from ..db_models.chado import Stock
from ..rest import Issue, IType, filter_parse, paginator


def create(**kwargs):
    content = None
    try:
        if not kwargs.get('uniquename'):
            raise Exception('Missing the uniquename')
        if not kwargs.get('type_id'):
            from biobarcoding.db_models.chado import Cvterm
            kwargs['type_id'] = chado_session.query(Cvterm.cvterm_id).filter(
                Cvterm.name == 'plant anatomical entity').one()
        chado_session.add(Stock(**kwargs))
        issues, status = [Issue(IType.INFO,
                                f'CREATE individuals: The individual "{kwargs.get("uniquename")}" created successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR,
                                f'CREATE individuals: The individual "{kwargs.get("uniquename")}" could not be created.')], 409
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
        issues, status = [Issue(IType.INFO, 'READ individuals: The individuals were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ individuals: The individuals could not be read.')], 400
    return issues, content, count, status


def update(id, **kwargs):
    content = None
    try:
        content = __get_query(id).first()
        content.update(kwargs)
        issues, status = [Issue(IType.INFO, f'UPDATE individuals: The individual "{id}" updated successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'UPDATE individuals: The individual "{id}" could not be updated.')], 409
    return issues, content, status


def delete(id=None, **kwargs):
    content = None
    try:
        query = __get_query(id, **kwargs)
        resp = query.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO,
                                f'DELETE individuals: The {resp} individuals were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE individuals: The individuals could not be removed.')], 404
    return issues, content, status


def __get_query(id=None, **kwargs):
    query = chado_session.query(Stock)
    global count
    count = 0
    if id:
        query = query.filter(Stock.stock_id == id)
    else:
        if 'filter' in kwargs:
            query = query.filter(filter_parse(Stock, kwargs.get('filter'), __aux_own_filter))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            count = query.count()
            query = paginator(query, kwargs.get('pagination'))
    return query


def __aux_own_filter(filter):
    clause = []

    if filter.get('feature_id'):
        from biobarcoding.db_models.chado import StockFeature
        _ids = chado_session.query(StockFeature.stock_id) \
            .filter(filter_parse(StockFeature, {'feature_id': filter.get('feature_id')}))
        clause.append(Stock.stock_id.in_(_ids))

    if 'phylotree_id' in filter:
        from biobarcoding.db_models.chado import Phylonode
        _ids = chado_session.query(Phylonode.feature_id) \
            .filter(
            filter_parse(Phylonode, [{'phylotree_id': filter.get('phylotree_id')}]))
        from biobarcoding.db_models.chado import StockFeature
        _ids = chado_session.query(StockFeature.stock_id) \
            .filter(filter_parse(StockFeature, {'feature_id': filter.get('feature_id')}))
        clause.append(Stock.stock_id.in_(_ids))

    if "cvterm_id" in filter:
        from biobarcoding.db_models.chado import StockCvterm
        _ids = chado_session.query(StockCvterm.stock_id) \
            .filter(filter_parse(StockCvterm, [{'cvterm_id': filter.get('cvterm_id')}]))
        clause.append(Stock.stock_id.in_(_ids))

    if "prop_cvterm_id" in filter:
        from biobarcoding.db_models.chado import Stockprop
        _ids = chado_session.query(Stockprop.stock_id) \
            .filter(filter_parse(Stockprop, [{'type_id': filter.get('prop_cvterm_id')}]))
        clause.append(Stock.stock_id.in_(_ids))

    # from datetime import datetime
    # if "added-from" in filter:
    #     filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
    #     _ids = chado_session.query(Stock.stock_id) \
    #         .filter(filter_parse(Stock, {'timeexecuted':filter.get("added-from")}))
    #     clause.append(Stock.stock_id.in_(_ids))
    # if "added-to" in filter:
    #     filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
    #     _ids = chado_session.query(Stock.stock_id) \
    #         .filter(filter_parse(Stock, {'timeexecuted':filter.get("added-to")}))
    #     clause.append(Stock.stock_id.in_(_ids))

    return clause


def __get_query_ordered(query, order):
    # query = query.order(order_parse(Stock, kwargs.get('order'), __aux_own_order))
    return query
