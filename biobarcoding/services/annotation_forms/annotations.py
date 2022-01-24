from .. import get_or_create
from ..main import SimpleAuxService, get_orm
from ...db_models import DBSession
from ...db_models.core import FunctionalObject
from ...db_models.sa_annotations import AnnotationItemFunctionalObject


##
# ANNOTATION TOOLS
##

class AuxService(SimpleAuxService):
    # TODO:
    #  if an existent field allow multiple instances
    #  if the value for a field is valid by its range
    #  if the value for a field is valid by its type (tag, attribute, relationship)

    def __init__(self):
        super(AuxService, self).__init__()
        self.orm = get_orm('annotations')

    def prepare_values(self, template=None, field=None, **values):

        if template:
            if not field and not values.get("type"):
                values["type"] = 'template'
            if not values.get('form_template') and not values.get('form_template_id'):
                from ...db_models.sa_annotations import AnnotationFormTemplate
                values['form_template'] = self.db.query(AnnotationFormTemplate)\
                    .filter(AnnotationFormTemplate.name == template).one()
        if field:
            if not template and not values.get("type"):
                values["type"] = 'field'
            if not values.get('form_field') and not values.get('form_field_id'):
                from ...db_models.sa_annotations import AnnotationFormField
                values['form_field'] = self.db.query(AnnotationFormField)\
                    .filter(AnnotationFormField.name == field).one()

        if values.get("type") in ('template', 'field', 'text'):
            self.orm = get_orm(f'annotation_{values.get("type")}')

        return super(AuxService, self).prepare_values(**values)

    def prepare_external_values(self, object_id=[], **values):
        if object_id and not values.get('object_uuid'):
            if isinstance(object_id, (tuple, list, set)):
                ids = object_id
            else:
                ids = [object_id]
            values['object_uuid'] = DBSession.query(FunctionalObject.uuid) \
                .filter(FunctionalObject.id.in_(ids)).all()
            # kwargs['object_uuid'] = [i for i, in kwargs['object_uuid']]
        return values

    def after_create(self, new_object, **values):
        if values.get('object_uuid'):
            if isinstance(values['object_uuid'], (tuple, list, set)):
                ids = values.get('object_uuid')
            else:
                ids = [values.get('object_uuid')]
            for i in ids:
                get_or_create(self.db, AnnotationItemFunctionalObject,
                              annotation_id=new_object.id, object_uuid=i)
        return values

    def after_update(self, new_object, **values):
        if values.get('object_uuid'):
            if isinstance(values['object_uuid'], (tuple, list, set)):
                ids = values.get('object_uuid')
            else:
                ids = [values.get('object_uuid')]
            rl = self.db.query(AnnotationItemFunctionalObject)\
                .filter(AnnotationItemFunctionalObject.annotation_id == new_object.id)
            rl.filter(AnnotationItemFunctionalObject.object_uuid.notin_(ids)) \
                .delete(synchronize_session='fetch')
            for i in ids:
                get_or_create(self.db, AnnotationItemFunctionalObject,
                              annotation_id=new_object.id, object_uuid=i)
        return values
