from . import MetaService
from ...main import get_orm
from ....db_models import DBSessionChado


##
# ANALYSIS SERVICE
##
class Service(MetaService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('analyses')

    ##
    # CREATE
    ##

    def prepare_values(self, **values):

        if values.get('job_id'):
            values['sourcename'] = values.get('job_id')
            values['sourceversion'] = 'job'
            values['sourceuri'] = f'/jobs/{values.get("job_id")}'

        return super(Service, self).prepare_values(**values)

    def check_values(self, **values):

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

        return super(Service, self).check_values(**values)

    ##
    # GETTER AND OTHERS
    ##

    def get_query(self, **kwargs):
        if kwargs.get('job_id'):
            # TODO: will there be multiple analyses for a single job ?
            query = self.db.query(self.orm).filter(self.orm.sourcename == str(kwargs.get('job_id')),
                                                   self.orm.sourceversion == 'job')
            return query, query.count()
        return super(Service, self).get_query(**kwargs)

    def aux_filter(self, filter):
        clauses = []
        from ....rest import filter_parse

        if filter.get('job_id'):
            _ids = self.db.query(self.orm.analysis_id) \
                .filter(filter_parse(self.orm, {'sourcename': filter.get('job_id'),
                                                'sourceversion': {'op': 'eq', 'unary': 'job'}}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get('feature_id'):
            from ....db_models.chado import AnalysisFeature
            _ids = self.db.query(AnalysisFeature.analysis_id) \
                .filter(filter_parse(AnalysisFeature, {'feature_id': filter.get('feature_id')}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get('organism_id'):
            from ....db_models.chado import Feature
            _ids = self.db.query(Feature.feature_id) \
                .filter(filter_parse(Feature, {'organism_id': filter.get('organism_id')}))
            from ....db_models.chado import AnalysisFeature
            _ids = self.db.query(AnalysisFeature.analysis_id) \
                .filter(AnalysisFeature.feature_id.in_(_ids))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if 'phylotree_id' in filter:
            from ....db_models.chado import Phylotree
            _ids = self.db.query(Phylotree.analysis_id) \
                .filter(filter_parse(Phylotree, [{'phylotree_id': filter.get('phylotree_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if "cvterm_id" in filter:
            from ....db_models.chado import AnalysisCvterm
            _ids = self.db.query(AnalysisCvterm.analysis_id) \
                .filter(filter_parse(AnalysisCvterm, [{'cvterm_id': filter.get('cvterm_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if "prop_cvterm_id" in filter:
            from ....db_models.chado import Analysisprop
            _ids = self.db.query(Analysisprop.analysis_id) \
                .filter(filter_parse(Analysisprop, [{'type_id': filter.get('prop_cvterm_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        from datetime import datetime
        if "added-from" in filter:
            filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
            _ids = self.db.query(self.orm.analysis_id) \
                .filter(filter_parse(self.orm, {'timeexecuted': filter.get("added-from")}))
            clauses.append(self.orm.analysis_id.in_(_ids))
        if "added-to" in filter:
            filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
            _ids = self.db.query(self.orm.analysis_id) \
                .filter(filter_parse(self.orm, {'timeexecuted': filter.get("added-to")}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(filter)
