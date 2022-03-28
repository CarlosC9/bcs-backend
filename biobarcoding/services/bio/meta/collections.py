from . import MetaService
from ...main import get_orm
from ....db_models import DBSessionChado


##
# COLLECTION SERVICE
##
class Service(MetaService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('collections')

    def check_values(self, **values):

        if not values.get('uniquename'):
            raise Exception('Missing the uniquename')

        if not values.get('type_id'):
            from ontologies import get_type_id
            values['type_id'] = get_type_id(type=values.get('type', 'stockcoll'))

        return super(Service, self).check_values(**values)

    def aux_filter(self, filter):
        from ....rest import filter_parse
        clauses = []

        if filter.get('stock_id'):
            from ....db_models.chado import StockcollectionStock
            _ids = self.db.query(StockcollectionStock.stockcollection_id) \
                .filter(filter_parse(StockcollectionStock, {'stock_id': filter.get('stock_id')}))
            clauses.append(self.orm.stockcollection_id.in_(_ids))

        if filter.get('feature_id'):
            from ....db_models.chado import Stock, StockcollectionStock
            _ids = self.db.query(Stock.stock_id) \
                .filter(filter_parse(Stock, {'feature_id': filter.get('feature_id')}))
            _ids = self.db.query(StockcollectionStock.stockcollection_id) \
                .filter(StockcollectionStock.stock_id.in_(_ids))
            clauses.append(self.orm.stockcollection_id.in_(_ids))

        if filter.get('organism_id'):
            from ....db_models.chado import Stock, StockcollectionStock
            _ids = self.db.query(Stock.stock_id) \
                .filter(filter_parse(Stock, {'organism_id': filter.get('organism_id')}))
            _ids = self.db.query(StockcollectionStock.stockcollection_id) \
                .filter(StockcollectionStock.stock_id.in_(_ids))
            clauses.append(self.orm.stockcollection_id.in_(_ids))

        if 'phylotree_id' in filter:
            # from ....db_models.chado import Phylotree
            # _ids = self.db.query(Phylotree.analysis_id) \
            #     .filter(filter_parse(Phylotree, [{'phylotree_id': filter.get('phylotree_id')}]))
            # clauses.append(self.orm.stockcollection_id.in_(_ids))
            pass

        if "cvterm_id" in filter:
            from ....db_models.chado import StockcollectionCvterm
            _ids = self.db.query(StockcollectionCvterm.stockcollection_id) \
                .filter(filter_parse(StockcollectionCvterm, [{'cvterm_id': filter.get('cvterm_id')}]))
            clauses.append(self.orm.stockcollection_id.in_(_ids))

        if "prop_cvterm_id" in filter:
            from ....db_models.chado import Stockcollectionprop
            _ids = self.db.query(Stockcollectionprop.stockcollection_id) \
                .filter(filter_parse(Stockcollectionprop, [{'type_id': filter.get('prop_cvterm_id')}]))
            clauses.append(self.orm.stockcollection_id.in_(_ids))

        # from datetime import datetime
        # if "added-from" in filter:
        #     filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
        #     _ids = self.db.query(self.orm.stockcollection_id) \
        #         .filter(filter_parse(self.orm, {'timeexecuted':filter.get("added-from")}))
        #     clauses.append(self.orm.stockcollection_id.in_(_ids))
        # if "added-to" in filter:
        #     filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
        #     _ids = self.db.query(self.orm.stockcollection_id) \
        #         .filter(filter_parse(self.orm, {'timeexecuted':filter.get("added-to")}))
        #     clauses.append(self.orm.stockcollection_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(filter)
