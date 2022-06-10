from .. import BioService
from ....rest import auth_filter, filter_parse
from ....db_models import DBSession
from ....db_models.core import data_object_type_id, FunctionalObject


##
# CHADO BOS SERVICE
##

class BosService(BioService):
    
    def __init__(self):
        super(BosService, self).__init__()
        self.obj_type = ''
        self.fos = FunctionalObject

    def prepare_values(self, **values):

        values['sourceuri'] = values.get('sourceuri', values.get('filesAPI'))

        return super(BosService, self).prepare_values(**values)

    def pre_query(self, purpose) -> object:
        from ....db_models.sysadmin import PermissionType
        purpose_id = DBSession.query(PermissionType).filter(PermissionType.name == purpose).one().id
        _ = data_object_type_id.get(self.obj_type)
        bos_clause = DBSession.query(self.fos.native_id) \
            .filter(auth_filter(self.fos, purpose_id, _ if not _ or isinstance(_, (tuple, list, set)) else [_]))
        bos_clause = [i for i, in bos_clause.all()]  # Cannot use .subquery(), because we have two -separate- databases

        from sqlalchemy import inspect
        bos_clause = inspect(self.orm).primary_key[0].in_(bos_clause)
        query = self.db.query(self.orm).filter(bos_clause)
        return query

    def aux_filter(self, filter):
        clauses = []

        if filter.get('annotation_field_id'):
            from ....db_models.sa_annotations import AnnotationItem, AnnotationItemFunctionalObject
            _ids = DBSession.query(self.fos.native_id) \
                .join(AnnotationItemFunctionalObject).join(AnnotationItem) \
                .filter(filter_parse(AnnotationItem, {'id': filter.get('annotation_field_id')}),
                        self.fos.obj_type_id == data_object_type_id[self.obj_type]).all()
            from sqlalchemy import inspect
            clauses.append(inspect(self.orm).primary_key[0].in_(_ids))

        if filter.get('annotation_form_template_id'):
            from ....db_models.sa_annotations import AnnotationFormTemplate, AnnotationTemplate, \
                AnnotationItemFunctionalObject
            _ids = DBSession.query(self.fos.native_id) \
                .join(AnnotationItemFunctionalObject).join(AnnotationTemplate).join(AnnotationFormTemplate) \
                .filter(filter_parse(AnnotationFormTemplate, {'id': filter.get('annotation_form_template_id')}),
                        self.fos.obj_type_id == data_object_type_id[self.obj_type]).all()
            from sqlalchemy import inspect
            clauses.append(inspect(self.orm).primary_key[0].in_(_ids))

        if filter.get('annotation_form_field_id'):
            from ....db_models.sa_annotations import AnnotationFormField, AnnotationField, \
                AnnotationItemFunctionalObject
            _ids = DBSession.query(self.fos.native_id) \
                .join(AnnotationItemFunctionalObject).join(AnnotationField).join(AnnotationFormField) \
                .filter(filter_parse(AnnotationFormField, {'id': filter.get('annotation_form_field_id')}),
                        self.fos.obj_type_id == data_object_type_id[self.obj_type]).all()
            from sqlalchemy import inspect
            clauses.append(inspect(self.orm).primary_key[0].in_(_ids))

        return clauses + super(BosService, self).aux_filter(filter)
