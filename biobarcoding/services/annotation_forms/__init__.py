from .. import get_or_create
from ..main import BasicService, get_orm
from ...rest import filter_parse
from ...db_models import DBSession, ObjectType
from ...db_models.sa_annotations import AnnotationFormItemObjectType


##
# ANNOTATION FORM SERVICE
##

class FormItemService(BasicService):

    def prepare_values(self, **values):

        if values.get("type"):
            self.orm = get_orm(f'{values.get("type")}') or self.orm

        return super(FormItemService, self).prepare_values(**values)

    def prepare_external_values(self, object_type=[], **values):
        if object_type and not values.get('object_type_id'):
            if isinstance(object_type, (tuple, list, set)):
                object_types = object_type
            else:
                object_types = [object_type]
            values['object_type_id'] = DBSession.query(ObjectType.id) \
                .filter(ObjectType.name.in_(object_types)).all()
            # kwargs['object_type_id'] = [i for i, in kwargs['object_type_id']]
        return super(FormItemService, self).prepare_external_values(**values)

    def after_create(self, new_object, **values):
        if values.get('object_type_id'):
            if isinstance(values['object_type_id'], (tuple, list, set)):
                ids = values.get('object_type_id')
            else:
                ids = [values.get('object_type_id')]
            for i in ids:
                get_or_create(self.db, AnnotationFormItemObjectType,
                              form_item_id=new_object.id, object_type_id=i)
        return values

    def read(self, **kwargs):
        content, count = self.get_query(purpose='read', **kwargs)
        if kwargs.get('id') or kwargs.get('cvterm_id') or kwargs.get('dbxref_id') \
                or (kwargs.get('cv') and kwargs.get('cvterm')) \
                or (kwargs.get('db') and kwargs.get('dbxref')):
            content = content.first()
        else:
            content = content.all()
        return content, count

    def after_update(self, new_object, **values):
        if values.get('object_type_id'):
            if isinstance(values['object_type_id'], (tuple, list, set)):
                ids = values.get('object_type_id')
            else:
                ids = [values.get('object_type_id')]
            rl = self.db.query(AnnotationFormItemObjectType) \
                .filter(AnnotationFormItemObjectType.form_item_id == new_object.id)
            rl.filter(AnnotationFormItemObjectType.object_type_id.notin_(ids)) \
                .delete(synchronize_session='fetch')
            for i in ids:
                get_or_create(self.db, AnnotationFormItemObjectType,
                              form_item_id=new_object.id, object_type_id=i)
        return values

    def aux_filter(self, _filter: dict) -> list:
        clauses = []

        if _filter.get('object_type') and not _filter.get('object_type_id'):
            _aux = DBSession.query(ObjectType.id).filter(
                filter_parse(ObjectType, {'name': _filter.get('object_type')}))
            _filter['object_type_id'] = {'op': 'in', 'unary': _aux}

        if _filter.get('object_type_id'):
            _aux = self.db.query(AnnotationFormItemObjectType.form_item_id).filter(
                filter_parse(AnnotationFormItemObjectType, {'object_type_id': _filter.get('object_type_id')}))
            clauses.append(self.orm.id.in_(_aux.subquery()))

        return clauses + super(FormItemService, self).aux_filter(_filter)