from .. import get_or_create, listify, get_filtering
from ..main import BasicService, get_orm
from ...db_models import DBSession
from ...db_models.core import FunctionalObject
from ...db_models.sa_annotations import AnnotationItemFunctionalObject


##
# ANNOTATION SERVICE
##

class Service(BasicService):
    # TODO:
    #  if an existent field allow multiple instances
    #  if the value for a field is valid by its range
    #  if the value for a field is valid by its type (tag, attribute, relationship)

    def __init__(self):
        super(Service, self).__init__()
        self.orm = get_orm('annotations')

    def prepare_values(self, template=None, field=None, **values):

        form_template = values.pop('form_template', None)
        template = form_template if not template and isinstance(form_template, str) else template
        if template:
            if not field and not values.get("type"):
                values["type"] = 'template'
            if not values.get('form_template') and not values.get('form_template_id'):
                from ...db_models.sa_annotations import AnnotationFormTemplate
                values['form_template'] = self.db.query(AnnotationFormTemplate) \
                    .filter(AnnotationFormTemplate.name == template).one()

        form_field = values.pop('form_field', None)
        field = form_field if not field and isinstance(form_field, str) else field
        if field:
            if not template and not values.get("type"):
                values["type"] = 'field'
            if not values.get('form_field') and not values.get('form_field_id'):
                from ...db_models.sa_annotations import AnnotationFormField
                values['form_field'] = self.db.query(AnnotationFormField) \
                    .filter(AnnotationFormField.name == field).one()

        if values.get("type") in ('template', 'field', 'text'):
            self.orm = get_orm(f'annotation_{values.get("type")}')

        return super(Service, self).prepare_values(**values)

    def prepare_external_values(self, object_id=None, **values):
        if object_id is not None and not values.get('object_uuid'):
            ids = listify(object_id)
            values['object_uuid'] = DBSession.query(FunctionalObject.uuid) \
                .filter(FunctionalObject.id.in_(ids)).all()
            # kwargs['object_uuid'] = [i for i, in kwargs['object_uuid']]
        return super(Service, self).prepare_external_values(**values)

    def after_create(self, new_object, **values):
        if values.get('object_uuid'):
            ids = listify(values.get('object_uuid'))
            for i in ids:
                rl = get_or_create(self.db, AnnotationItemFunctionalObject,
                                   annotation_id=new_object.id, object_uuid=i)
                if isinstance(values.get('rank'), int):
                    rl.rank = values.get('rank')
        return values

    def update(self, id=None, object_uuid=None, values={}, **kwargs):
        if id and not isinstance(values, (list, tuple)):
            return super(Service, self).update(id=id, values=values, **kwargs)
        elif object_uuid and isinstance(values, (list, tuple)):
            # DELETE anything missing
            self.db.query(AnnotationItemFunctionalObject) \
                .filter(AnnotationItemFunctionalObject.object_uuid == object_uuid) \
                .delete(synchronize_session='fetch')
            # UPDATE anything with id
            # CREATE anything without id
            content, count, i = [], 0, 0
            for v in values:
                id = v.get('id')
                i+=1
                if id:
                    c, cc = super(Service, self).update(id=id, values=v, **kwargs)
                    if len(c) != 1:
                        raise Exception(f'UPDATE: Could not be updated the annotation {id}')
                    self.after_create(c[0], object_uuid=object_uuid, rank=v.get('rank', i))
                    content.append(c)
                    count += cc
                else:
                    v['object_uuid'] = v.get('object_uuid', object_uuid)
                    c, cc = super(Service, self).create(**v)
                    content.append(c)
                    count += cc
            return content, count
        else:
            raise Exception('UPDATE: Bad request for update annotations.')

    def after_update(self, new_object, **values) -> dict:

        if values.get('object_uuid'):
            ids = listify(values.get('object_uuid'))
            rl = self.db.query(AnnotationItemFunctionalObject) \
                .filter(AnnotationItemFunctionalObject.annotation_id == new_object.id)
            rl.filter(AnnotationItemFunctionalObject.object_uuid.notin_(ids)) \
                .delete(synchronize_session='fetch')
            for i in ids:
                get_or_create(self.db, AnnotationItemFunctionalObject,
                              annotation_id=new_object.id, object_uuid=i)

        if values.get('new_object_uuid'):
            ids = listify(values.get('new_object_uuid'))
            for i in ids:
                get_or_create(self.db, AnnotationItemFunctionalObject,
                              annotation_id=new_object.id, object_uuid=i)

        return values

    def get_query(self, query=None, purpose='delete', object_uuid=None, **kwargs) -> (object, int):

        if get_filtering("type", kwargs) in ('template', 'field', 'text'):
            self.orm = get_orm(f'annotation_{get_filtering("type", kwargs)}')

        if object_uuid:
            query = query or self.pre_query(purpose) or self.db.query(self.orm)
            query = query.join(AnnotationItemFunctionalObject).filter(
                AnnotationItemFunctionalObject.object_uuid == object_uuid).order_by(AnnotationItemFunctionalObject.rank)

        return super(Service, self).get_query(query=query, **kwargs)

    def aux_filter(self, _filter: dict) -> list:
        clauses = []

        if _filter.get('object_uuid'):
            from ...rest import filter_parse
            _aux = self.db.query(AnnotationItemFunctionalObject.annotation_id).filter(
                filter_parse(AnnotationItemFunctionalObject, {'object_uuid': _filter.get('object_uuid')}))
            clauses.append(self.orm.id.in_(_aux.subquery()))

        return clauses + super(Service, self).aux_filter(_filter)
