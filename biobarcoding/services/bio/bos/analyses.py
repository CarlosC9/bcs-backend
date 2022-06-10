from . import BosService
from ..meta.ontologies import get_type_id
from ... import get_or_create
from ...main import get_orm
from ....db_models import DBSessionChado, DBSession
from ....db_models.chado import AnalysisRelationship


##
# ANALYSIS SERVICE
##
class Service(BosService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('analyses')
        self.obj_type = 'analysis'

    ##
    # CREATE
    ##

    def prepare_values(self, **values):

        if values.get('job_id'):
            from ....db_models.jobs import Job
            job = DBSession.query(Job).filter(Job.id == values.get('job_id')).one()
            # DBSession:Process
            values['programname'] = job.process_id
            values['programversion'] = job.process.name
            # DBSession:Job
            values['sourcename'] = job.id
            values['sourceversion'] = 'job'
            values['sourceuri'] = f'/jobs/{job.id}'

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

    def after_create(self, new_object, **values):
        values = super(Service, self).after_create(new_object, **values)

        fos = get_or_create(DBSession, self.fos,
                            native_id=new_object.analysis_id,
                            name=new_object.name)

        if values.get('derives_from'):
            _ = values.get('derives_from')
            _ = _ if isinstance(_, (tuple, list, set)) else [_]
            rls = [(i, new_object.analysis_id) for i in _]
        elif values.get('derives_into'):
            _ = values.get('derives_into')
            _ = _ if isinstance(_, (tuple, list, set)) else [_]
            rls = [(new_object.analysis_id, i) for i in _]
        else:
            return values
        for sbj, obj in rls:
            get_or_create(self.db, AnalysisRelationship,
                          subject_id=sbj, object_id=obj,
                          type_id=get_type_id(type='relationship', subtype='derives_into'))

        return values

    ##
    # DELETE
    ##

    def delete_related(self, *content, **kwargs):
        # TODO: delete from Phylotree ? check cascade with Analysis
        ids = [t.analysis_id for t in content]
        query = DBSession.query(self.fos).filter(self.fos.native_id.in_(ids))
        return len([DBSession.delete(row) for row in query.all()])

    ##
    # GETTER AND OTHERS
    ##

    def aux_filter(self, filter):
        clauses = []
        from ....rest import filter_parse

        if filter.get('job_id'):
            _ids = self.db.query(self.orm.analysis_id) \
                .filter(filter_parse(self.orm, {'sourcename': filter.get('job_id'),
                                                'sourceversion': 'job'}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get('derives_from'):
            _ids = self.db.query(AnalysisRelationship.object_id) \
                .filter(filter_parse(AnalysisRelationship,
                                     {'subject_id': filter.get('derives_from'),
                                      'type_id': get_type_id(type='relationship', subtype='derives_into')}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get('derives_into'):
            _ids = self.db.query(AnalysisRelationship.subject_id) \
                .filter(filter_parse(AnalysisRelationship,
                                     {'object_id': filter.get('derives_into'),
                                      'type_id': get_type_id(type='relationship', subtype='derives_into')}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get('feature_id'):
            from ....db_models.chado import AnalysisFeature
            _ids = self.db.query(AnalysisFeature.analysis_id) \
                .filter(filter_parse(AnalysisFeature, {'feature_id': filter.get('feature_id')}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get('organism_id'):
            from ....db_models.chado import AnalysisFeature, Feature
            _ids = self.db.query(AnalysisFeature.analysis_id).join(Feature) \
                .filter(filter_parse(Feature, {'organism_id': filter.get('organism_id')}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get('phylotree_id'):
            from ....db_models.chado import Phylotree
            _ids = self.db.query(Phylotree.analysis_id) \
                .filter(filter_parse(Phylotree, [{'phylotree_id': filter.get('phylotree_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get("cvterm_id"):
            from ....db_models.chado import AnalysisCvterm
            _ids = self.db.query(AnalysisCvterm.analysis_id) \
                .filter(filter_parse(AnalysisCvterm, [{'cvterm_id': filter.get('cvterm_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get("prop_cvterm_id"):
            from ....db_models.chado import Analysisprop
            _ids = self.db.query(Analysisprop.analysis_id) \
                .filter(filter_parse(Analysisprop, [{'type_id': filter.get('prop_cvterm_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        from datetime import datetime
        if filter.get("added-from"):
            filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
            _ids = self.db.query(self.orm.analysis_id) \
                .filter(filter_parse(self.orm, {'timeexecuted': filter.get("added-from")}))
            clauses.append(self.orm.analysis_id.in_(_ids))
        if filter.get("added-to"):
            filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
            _ids = self.db.query(self.orm.analysis_id) \
                .filter(filter_parse(self.orm, {'timeexecuted': filter.get("added-to")}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(filter)
