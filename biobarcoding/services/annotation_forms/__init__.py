from .. import get_or_create, listify
from ..main import BasicService, get_orm
from ...rest import filter_parse
from ...db_models import DBSession, ObjectType, DBSessionChado
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
            object_types = listify(object_type)
            values['object_type_id'] = DBSession.query(ObjectType.id) \
                .filter(ObjectType.name.in_(object_types)).all()
            # kwargs['object_type_id'] = [i for i, in kwargs['object_type_id']]
        return super(FormItemService, self).prepare_external_values(**values)

    def after_create(self, new_object, **values):
        if values.get('object_type_id'):
            ids = listify(values.get('object_type_id'))
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
            ids = listify(values.get('object_type_id'))
            rl = self.db.query(AnnotationFormItemObjectType) \
                .filter(AnnotationFormItemObjectType.form_item_id == new_object.id)
            rl.filter(AnnotationFormItemObjectType.object_type_id.notin_(ids)) \
                .delete(synchronize_session='fetch')
            for i in ids:
                get_or_create(self.db, AnnotationFormItemObjectType,
                              form_item_id=new_object.id, object_type_id=i)
        return values

    def get_query(self, query=None, id=None, purpose='delete', **kwargs) -> (object, int):
        if kwargs.get('db') and kwargs.get('dbxref'):
            from ...db_models.chado import Db, Dbxref
            kwargs['dbxref_id'] = kwargs.get('dbxref_id') or DBSessionChado.query(Dbxref.dbxref_id).join(Db) \
                .filter(Db.name == kwargs.get('db'), Dbxref.accession == kwargs.get('dbxref')).one()
        return super(FormItemService, self).get_query(query, id, purpose, **kwargs)

    def aux_filter(self, _filter: dict) -> list:
        clauses = []

        if _filter.get('object_type') and not _filter.get('object_type_id'):
            _aux = DBSession.query(ObjectType.id).filter(
                filter_parse(ObjectType, {'name': _filter.get('object_type')}))
            _filter['object_type_id'] = {'op': 'in', 'unary': _aux.subquery()}

        if _filter.get('object_type_id'):
            _aux = self.db.query(AnnotationFormItemObjectType.form_item_id).filter(
                filter_parse(AnnotationFormItemObjectType, {'object_type_id': _filter.get('object_type_id')}))
            clauses.append(self.orm.id.in_(_aux.subquery()))

        return clauses + super(FormItemService, self).aux_filter(_filter)