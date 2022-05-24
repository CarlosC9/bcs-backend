from . import MetaService
from ...main import get_orm
from ....db_models import DBSession, DBSessionChado
from ....db_models.bioinformatics import Specimen
from ... import get_or_create


##
# INDIVIDUAL SERVICE
##
class Service(MetaService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('individuals')

    ##
    # CREATE
    ##

    def check_values(self, **values):

        if not values.get('uniquename'):
            raise Exception('Missing the uniquename')

        if not values.get('type_id'):
            from .ontologies import get_type_id
            values['type_id'] = get_type_id(type=values.get('type', 'stock'))

        return super(Service, self).check_values(**values)

    def after_create(self, stock, **values):
        super(Service, self).after_create(stock, **values)

        s = get_or_create(DBSession, Specimen,
                          # native_id=stock.stock_id,
                          # native_table='stock',
                          name=stock.uniquename)

        return values

    ##
    # DELETE
    ##

    def delete_related(self, *content, **kwargs):
        names = [s.uniquename for s in content]
        query = DBSession.query(Specimen).filter(Specimen.name.in_(names))
        return len([DBSession.delete(row) for row in query.all()])

    ##
    # GET SQLALCHEMY QUERY
    ##

    def aux_filter(self, filter):
        from ....rest import filter_parse
        clauses = []

        if filter.get('feature_id'):
            from ....db_models.chado import StockFeature
            _ids = self.db.query(StockFeature.stock_id) \
                .filter(filter_parse(StockFeature, {'feature_id': filter.get('feature_id')}))
            clauses.append(self.orm.stock_id.in_(_ids))

        if 'phylotree_id' in filter:
            from ....db_models.chado import Phylonode
            _ids = self.db.query(Phylonode.feature_id) \
                .filter(
                filter_parse(Phylonode, [{'phylotree_id': filter.get('phylotree_id')}]))
            from ....db_models.chado import StockFeature
            _ids = self.db.query(StockFeature.stock_id) \
                .filter(filter_parse(StockFeature, {'feature_id': filter.get('feature_id')}))
            clauses.append(self.orm.stock_id.in_(_ids))

        if "cvterm_id" in filter:
            from ....db_models.chado import StockCvterm
            _ids = self.db.query(StockCvterm.stock_id) \
                .filter(filter_parse(StockCvterm, [{'cvterm_id': filter.get('cvterm_id')}]))
            clauses.append(self.orm.stock_id.in_(_ids))

        if "prop_cvterm_id" in filter:
            from ....db_models.chado import Stockprop
            _ids = self.db.query(Stockprop.stock_id) \
                .filter(filter_parse(Stockprop, [{'type_id': filter.get('prop_cvterm_id')}]))
            clauses.append(self.orm.stock_id.in_(_ids))

        # from datetime import datetime
        # if "added-from" in filter:
        #     filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
        #     _ids = self.db.query(self.orm.stock_id) \
        #         .filter(filter_parse(self.orm, {'timeexecuted':filter.get("added-from")}))
        #     clauses.append(self.orm.stock_id.in_(_ids))
        # if "added-to" in filter:
        #     filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
        #     _ids = self.db.query(self.orm.stock_id) \
        #         .filter(filter_parse(self.orm, {'timeexecuted':filter.get("added-to")}))
        #     clauses.append(self.orm.stock_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(filter)
