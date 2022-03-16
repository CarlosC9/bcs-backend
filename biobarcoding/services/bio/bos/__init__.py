from ... import get_orm_params, get_or_create
from ...main import SimpleAuxService
from ....rest import filter_parse
from ....db_models import DBSession, DBSessionChado, ObjectType
from ....db_models.sysadmin import Collection
from ....db_models.chado import Cv, Cvterm, Db, Dbxref
from ....db_models.sa_annotations import AnnotationFormItemObjectType


##
# BOS TOOLS
##

class BosAuxService(SimpleAuxService):

    # def prepare_values(self, cv=None, cvterm=None, db=None, dbxref=None, **values):
    #
    #     if cv and cvterm and not values.get('cvterm_id'):
    #         values['cvterm_id'] = DBSessionChado.query(Cvterm).join(Cv) \
    #             .filter(Cv.name == cv, Cvterm.name == cvterm).one().cvterm_id
    #
    #     if db and dbxref and not values.get('dbxref_id'):
    #         values['dbxref_id'] = DBSessionChado.query(Dbxref).join(Db) \
    #             .filter(Db.name == db, Dbxref.accession == dbxref).one().dbxref_id
    #
    #     if values.get('cvterm_id') and not values.get('dbxref_id'):
    #         values['dbxref_id'] = DBSessionChado.query(Cvterm) \
    #             .filter(Cvterm.cvterm_id == values.get('cvterm_id'))\
    #             .first().dbxref_id
    #
    #     if values.get('dbxref_id') and not values.get('cvterm_id'):
    #         values['cvterm_id'] = DBSessionChado.query(Cvterm) \
    #             .filter(Cvterm.dbxref_id == values.get('dbxref_id'))\
    #             .first().cvterm_id
    #
    #     return get_orm_params(self.orm, **values)
    #
    # def prepare_external_values(self, object_type=[], **values):
    #
    #     if object_type and not values.get('object_type_id'):
    #         if isinstance(object_type, (tuple, list, set)):
    #             object_types = object_type
    #         else:
    #             object_types = [object_type]
    #         values['object_type_id'] = DBSession.query(ObjectType.object_type_id) \
    #             .filter(ObjectType.name.in_(object_types)).all()
    #         # kwargs['object_type_id'] = [i for i, in kwargs['object_type_id']]
    #
    #     return values
    #
    # def after_create(self, new_object, object_type=[], **kwargs):
    #     values = self.prepare_external_values(**kwargs)
    #     if values.get('object_type_id'):
    #         if isinstance(values['object_type_id'], (tuple, list, set)):
    #             ids = values.get('object_type_id')
    #         else:
    #             ids = [values.get('object_type_id')]
    #         for i in ids:
    #             get_or_create(DBSession, AnnotationFormItemObjectType,
    #                           form_item_id=new_object.id, object_type_id=i)
    #     return None

    def read(self, key_id, **kwargs):
        content, count = self.get_query(**kwargs)
        if kwargs.get(key_id):
            content = content.first()
        else:
            content = content.all()
        return content, count

    def aux_filter(self, filter):
        clauses = []

        if filter.get('collection') and not filter.get('collection_id'):
            filter['collection_id'] = DBSession.query(Collection.collection_id) \
                .filter(Collection.name.in_(filter.get('collection'))).all()
            filter['collection_id'] = {'op':'in', 'unary': [i for i, in filter['collection_id']]}

        if filter.get('collection_id'):
            _ids = DBSession.query(AnnotationFormItemObjectType.form_item_id) \
                .filter(filter_parse(AnnotationFormItemObjectType, {'collection_id': filter.get('collection_id')}))
            clauses.append(self.orm.id.in_(_ids.subquery()))

        return clauses