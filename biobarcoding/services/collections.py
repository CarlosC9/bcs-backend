from ..db_models import DBSessionChado as chado_session
from ..db_models.chado import Stockcollection
from ..rest import Issue, IType, filter_parse
from . import get_query


def create(**kwargs):
    content = None
    try:
        if not kwargs.get('uniquename'):
            raise Exception('Missing the uniquename')
        if not kwargs.get('type_id'):
            from biobarcoding.db_models.chado import Cvterm
            kwargs['type_id'] = chado_session.query(Cvterm.cvterm_id).filter(Cvterm.name == 'sequence_collection').one()
        chado_session.add(Stockcollection(**kwargs))
        issues, status = [Issue(IType.INFO,
                                f'CREATE collections: The collection "{kwargs.get("uniquename")}" created successfully.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR,
                                f'CREATE collections: The collection "{kwargs.get("uniquename")}" could not be created.')], 409
    return issues, content, status


def read(id=None, **kwargs):
    content, count = None, 0
    try:
        content, count = __get_query(id, **kwargs)
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
        content, count = __get_query(id)
        content = content.update(kwargs)
        issues, status = [Issue(IType.INFO, f'UPDATE collections: The collection "{id}" was successfully updated.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'UPDATE collections: The collection "{id}" could not be updated.')], 409
    return issues, content, status


def delete(id=None, **kwargs):
    content = None
    try:
        content, count = __get_query(id, **kwargs)
        content = content.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO,
                                f'DELETE collections: The {content} collections were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE collections: The collections could not be removed.')], 404
    return issues, content, status


def __get_query(stockcollection_id=None, **kwargs):
    if stockcollection_id:
        query = chado_session.query(Stockcollection).filter(Stockcollection.stockcollection_id == stockcollection_id)
        return query, query.count()
    return get_query(chado_session, Stockcollection, **kwargs,
                     aux_filter=__aux_own_filter, aux_order=__aux_own_order)


def __aux_own_filter(filter):
    clause = []

    if filter.get('stock_id'):
        from biobarcoding.db_models.chado import StockcollectionStock
        _ids = chado_session.query(StockcollectionStock.stockcollection_id) \
            .filter(filter_parse(StockcollectionStock, {'stock_id': filter.get('stock_id')}))
        clause.append(Stockcollection.stockcollection_id.in_(_ids))

    if filter.get('feature_id'):
        from biobarcoding.db_models.chado import Stock, StockcollectionStock
        _ids = chado_session.query(Stock.stock_id) \
            .filter(filter_parse(Stock, {'feature_id': filter.get('feature_id')}))
        _ids = chado_session.query(StockcollectionStock.stockcollection_id) \
            .filter(StockcollectionStock.stock_id.in_(_ids))
        clause.append(Stockcollection.stockcollection_id.in_(_ids))

    if filter.get('organism_id'):
        from biobarcoding.db_models.chado import Stock, StockcollectionStock
        _ids = chado_session.query(Stock.stock_id) \
            .filter(filter_parse(Stock, {'organism_id': filter.get('organism_id')}))
        _ids = chado_session.query(StockcollectionStock.stockcollection_id) \
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


def __aux_own_order(order):
    # query = query.order(order_parse(Stockcollection, kwargs.get('order'), __aux_own_order))
    return []
