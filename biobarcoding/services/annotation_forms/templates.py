from . import FormItemAuxService
from ..main import get_orm
from ...rest import filter_parse
from ...services import get_or_create
from ...db_models.sa_annotations import AnnotationFormField, AnnotationFormTemplateField


##
# TEMPLATE TOOLS
##

class AuxService(FormItemAuxService):

    def __init__(self):
        super(AuxService, self).__init__()
        self.orm = get_orm('templates')

    def prepare_external_values(self, field=None, required_field=None, **values):

        if field is not None and not values.get('field_id'):
            if isinstance(field, (tuple, list, set)):
                filter = AnnotationFormField.name.in_(field)
            else:
                filter = AnnotationFormField.name == field
            values['field_id'] = self.db.query(AnnotationFormField.id) \
                .filter(filter).all()

        if required_field is not None and not values.get('required_field_id'):
            if isinstance(required_field, (tuple, list, set)):
                filter = AnnotationFormField.name.in_(required_field)
            else:
                filter = AnnotationFormField.name == required_field
            values['required_field_id'] = self.db.query(AnnotationFormField.id) \
                .filter(filter).all()

        return super(AuxService, self).prepare_external_values(**values)

    def after_create(self, new_object, **values):
        super(AuxService, self).after_create(new_object, **values)

        if values.get('field_id'):
            if isinstance(values['field_id'], (tuple, list, set)):
                ids = values.get('field_id')
            else:
                ids = [values.get('field_id')]
            for i in ids:
                get_or_create(self.db, AnnotationFormTemplateField,
                              form_template=new_object, form_field_id=i)

        if values.get('required_field_id'):
            if isinstance(values['required_field_id'], (tuple, list, set)):
                ids = values.get('required_field_id')
            else:
                ids = [values.get('required_field_id')]
            for i in ids:
                get_or_create(self.db, AnnotationFormTemplateField,
                              form_template=new_object, form_field_id=i, required=True)

        return values

    def after_update(self, new_object, **values):
        super(AuxService, self).after_update(new_object, **values)

        if values.get('field_id') is not None:
            if isinstance(values['field_id'], (tuple, list, set)):
                ids = values.get('field_id')
            else:
                ids = [values.get('field_id')]
            rl = self.db.query(AnnotationFormTemplateField) \
                .filter(AnnotationFormTemplateField.form_template_id == new_object.id,
                        AnnotationFormTemplateField.required is False)
            rl.filter(AnnotationFormTemplateField.form_field_id.notin_(ids)) \
                .delete(synchronize_session='fetch')
            for i in ids:
                get_or_create(self.db, AnnotationFormTemplateField,
                              form_template_id=new_object.id, form_field_id=i)

        if values.get('required_field_id') is not None:
            if isinstance(values['required_field_id'], (tuple, list, set)):
                ids = values.get('required_field_id')
            else:
                ids = [values.get('required_field_id')]
            rl = self.db.query(AnnotationFormTemplateField) \
                .filter(AnnotationFormTemplateField.form_template_id == new_object.id,
                        AnnotationFormTemplateField.required is True)
            rl.filter(AnnotationFormTemplateField.form_field_id.notin_(ids)) \
                .delete(synchronize_session='fetch')
            for i in ids:
                get_or_create(self.db, AnnotationFormTemplateField,
                              form_template_id=new_object.id, form_field_id=i, required=True)

        return values

    def aux_filter(self, filter):
        clauses = []

        if filter.get('field_id'):
            _ids = self.db.query(AnnotationFormTemplateField.form_template_id) \
                .filter(filter_parse(AnnotationFormTemplateField, {'field_id': filter.get('field_id')}))
            clauses.append(self.orm.id.in_(_ids.subquery()))

        return clauses + super(AuxService, self).aux_filter(filter)