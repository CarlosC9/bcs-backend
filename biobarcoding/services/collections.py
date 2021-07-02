from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.db_models.chado import Stockcollection
from biobarcoding.rest import Issue, IType, filter_parse, paginator


def create(**kwargs):
    content = None
    try:
        if not kwargs.get('uniquename'):
            raise Exception('Missing the uniquename')
        if not kwargs.get('type_id'):
            from biobarcoding.db_models.chado import Cvterm
            kwargs['type_id'] = chado_session.query(Cvterm.cvterm_id).filter(Cvterm.name=='sequence_collection').one()
        chado_session.add(Stockcollection(**kwargs))
        issues, status = [Issue(IType.INFO, f'CREATE collections: The collection "{kwargs.get("uniquename")}" created successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE collections: The collection "{kwargs.get("uniquename")}" could not be created.')], 409
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
        issues, status = [Issue(IType.INFO, 'READ collections: The collections were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ collections: The collections could not be read.')], 400
    return issues, content, count, status


def update(id, **kwargs):
    content = None
    try:
        content = __get_query(id).update(kwargs)
        issues, status = [Issue(IType.INFO, f'UPDATE collections: The collection "{id}" updated successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'UPDATE collections: The collection "{id}" could not be updated.')], 409
    return issues, content, status


def delete(id=None, **kwargs):
    content = None
    try:
        query = __get_query(id, **kwargs)
        resp = query.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE collections: The {resp} collections were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE collections: The collections could not be removed.')], 404
    return issues, content, status


def __get_query(id=None, **kwargs):
    query = chado_session.query(Stockcollection)
    global count
    count = 0
    if id:
        query = query.filter(Stockcollection.stockcollection_id == id)
    else:
        if 'filter' in kwargs:
            query = query.filter(filter_parse(Stockcollection, kwargs.get('filter'), __aux_own_filter))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            count = query.count()
            query = paginator(query, kwargs.get('pagination'))
    return query


def __aux_own_filter(filter):
    clause = []

    if filter.get('stock_id'):
        from biobarcoding.db_models.chado import StockcollectionStock
        _ids = chado_session.query(StockcollectionStock.stockcollection_id)\
            .filter(filter_parse(StockcollectionStock, {'stock_id': filter.get('stock_id')}))
        clause.append(Stockcollection.stockcollection_id.in_(_ids))

    if filter.get('feature_id'):
        from biobarcoding.db_models.chado import Stock, StockcollectionStock
        _ids = chado_session.query(Stock.stock_id)\
            .filter(filter_parse(Stock, {'feature_id': filter.get('feature_id')}))
        _ids = chado_session.query(StockcollectionStock.stockcollection_id)\
            .filter(StockcollectionStock.stock_id.in_(_ids))
        clause.append(Stockcollection.stockcollection_id.in_(_ids))

    if filter.get('organism_id'):
        from biobarcoding.db_models.chado import Stock, StockcollectionStock
        _ids = chado_session.query(Stock.stock_id)\
            .filter(filter_parse(Stock, {'organism_id': filter.get('organism_id')}))
        _ids = chado_session.query(StockcollectionStock.stockcollection_id)\
            .filter(StockcollectionStock.stock_id.in_(_ids))
        clause.append(Stockcollection.stockcollection_id.in_(_ids))

    if 'phylotree_id' in filter:
        # from biobarcoding.db_models.chado import Phylotree
        # _ids = chado_session.query(Phylotree.analysis_id) \
        #     .filter(filter_parse(Phylotree, [{'phylotree_id': filter.get('phylotree_id')}]))
        # clause.append(Stockcollection.stockcollection_id.in_(_ids))
        pass

    if "cvterm_id" in filter:
        from biobarcoding.db_models.chado import StockcollectionCvterm
        _ids = chado_session.query(StockcollectionCvterm.stockcollection_id) \
            .filter(filter_parse(StockcollectionCvterm, [{'cvterm_id': filter.get('cvterm_id')}]))
        clause.append(Stockcollection.stockcollection_id.in_(_ids))

    if "prop_cvterm_id" in filter:
        from biobarcoding.db_models.chado import Stockcollectionprop
        _ids = chado_session.query(Stockcollectionprop.stockcollection_id) \
            .filter(filter_parse(Stockcollectionprop, [{'type_id': filter.get('prop_cvterm_id')}]))
        clause.append(Stockcollection.stockcollection_id.in_(_ids))

    # from datetime import datetime
    # if "added-from" in filter:
    #     filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
    #     _ids = chado_session.query(Stockcollection.stockcollection_id) \
    #         .filter(filter_parse(Stockcollection, {'timeexecuted':filter.get("added-from")}))
    #     clause.append(Stockcollection.stockcollection_id.in_(_ids))
    # if "added-to" in filter:
    #     filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
    #     _ids = chado_session.query(Stockcollection.stockcollection_id) \
    #         .filter(filter_parse(Stockcollection, {'timeexecuted':filter.get("added-to")}))
    #     clause.append(Stockcollection.stockcollection_id.in_(_ids))

    return clause


def __get_query_ordered(query, order):
    # query = query.order(order_parse(Stockcollection, kwargs.get('order'), __aux_own_order))
    return query