from .. import get_or_create
from ..main import BasicService, get_orm
from ...rest import filter_parse
from ...db_models import DBSession, DBSessionChado, ObjectType
from ...db_models.chado import Cv, Cvterm, Db, Dbxref
from ...db_models.sa_annotations import AnnotationFormItemObjectType


##
# ANNOTATION FORM SERVICE
##

class FormItemService(BasicService):

    def prepare_values(self, cv=None, cvterm=None, db=None, dbxref=None, **values):
        if cv and cvterm and not values.get('cvterm_id'):
            from sqlalchemy.orm.exc import NoResultFound
            try:
                values['cvterm_id'] = DBSessionChado.query(Cvterm).join(Cv) \
                    .filter(Cv.name == cv, Cvterm.name == cvterm).one().cvterm_id
            except NoResultFound as e:
                pass

        if db and dbxref and not values.get('dbxref_id'):
            values['dbxref_id'] = DBSessionChado.query(Dbxref).join(Db) \
                .filter(Db.name == db, Dbxref.accession == dbxref).one().dbxref_id

        if values.get('cvterm_id') and not values.get('dbxref_id'):
            values['dbxref_id'] = DBSessionChado.query(Cvterm) \
                .filter(Cvterm.cvterm_id == values.get('cvterm_id')) \
                .first().dbxref_id

        if values.get('dbxref_id') and not values.get('cvterm_id'):
            values['cvterm_id'] = DBSessionChado.query(Cvterm) \
                .filter(Cvterm.dbxref_id == values.get('dbxref_id')) \
                .first().cvterm_id

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