from . import MetaService
from ...main import get_orm
from ....db_models import DBSession


##
# COLLECTION SERVICE
##
class Service(MetaService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSession
        self.orm = get_orm('collections')


##
# CHADO STOCK COLLECTION SERVICE
##
class StockCollService(MetaService):

    def __init__(self):
        super(StockCollService, self).__init__()
        from ....db_models.chado import Stockcollection
        self.orm = Stockcollection

    def check_values(self, **values):

        if not values.get('uniquename'):
            raise Exception('Missing the uniquename')

        if not values.get('type_id'):
            from .ontologies import get_type_id
            values['type_id'] = get_type_id(type=values.get('type', 'stockcoll'))

        return super(StockCollService, self).check_values(**values)

    def aux_filter(self, _filter: dict) -> list:
        from ....rest import filter_parse
        clauses = []

        if _filter.get('stock_id'):
            from ....db_models.chado import StockcollectionStock
            _ids = self.db.query(StockcollectionStock.stockcollection_id) \
                .filter(filter_parse(StockcollectionStock, {'stock_id': _filter.get('stock_id')}))
            clauses.append(self.orm.stockcollection_id.in_(_ids))

        if _filter.get('feature_id'):
            from ....db_models.chado import Stock, StockcollectionStock
            _ids = self.db.query(StockcollectionStock.stockcollection_id).join(Stock) \
                .filter(filter_parse(Stock, {'feature_id': _filter.get('feature_id')}))
            clauses.append(self.orm.stockcollection_id.in_(_ids))

        if _filter.get('organism_id'):
            from ....db_models.chado import Stock, StockcollectionStock
            _ids = self.db.query(StockcollectionStock.stockcollection_id).join(Stock) \
                .filter(filter_parse(Stock, {'organism_id': _filter.get('organism_id')}))
            clauses.append(self.orm.stockcollection_id.in_(_ids))

        if _filter.get("cvterm_id"):
            from ....db_models.chado import StockcollectionCvterm
            _ids = self.db.query(StockcollectionCvterm.stockcollection_id) \
                .filter(filter_parse(StockcollectionCvterm, [{'cvterm_id': _filter.get('cvterm_id')}]))
            clauses.append(self.orm.stockcollection_id.in_(_ids))

        if _filter.get("prop_cvterm_id"):
            from ....db_models.chado import Stockcollectionprop
            _ids = self.db.query(Stockcollectionprop.stockcollection_id) \
                .filter(filter_parse(Stockcollectionprop, [{'type_id': _filter.get('prop_cvterm_id')}]))
            clauses.append(self.orm.stockcollection_id.in_(_ids))

        # from datetime import datetime
        # if filter.get("added-from"):
        #     filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
        #     _ids = self.db.query(self.orm.stockcollection_id) \
        #         .filter(filter_parse(self.orm, {'timeexecuted':filter.get("added-from")}))
        #     clauses.append(self.orm.stockcollection_id.in_(_ids))
        # if filter.get("added-to"):
        #     filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
        #     _ids = self.db.query(self.orm.stockcollection_id) \
        #         .filter(filter_parse(self.orm, {'timeexecuted':filter.get("added-to")}))
        #     clauses.append(self.orm.stockcollection_id.in_(_ids))

        return clauses + super(StockCollService, self).aux_filter(_filter)
