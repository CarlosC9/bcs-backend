from . import FormItemAuxService
from ..main import get_orm
from ...rest import filter_parse
from ...services import get_or_create
from ...db_models.sa_annotations import AnnotationFormTemplate, AnnotationFormTemplateField


##
# FIELD TOOLS
##

class AuxService(FormItemAuxService):

    def __init__(self):
        super(AuxService, self).__init__()
        self.orm = get_orm('fields')

    def prepare_external_values(self, template=[], **values):
        if template and not values.get('template_id'):
            if isinstance(template, (tuple, list, set)):
                filter = AnnotationFormTemplate.name.in_(template)
            else:
                filter = AnnotationFormTemplate.name == template
            values['template_id'] = self.db.query(AnnotationFormTemplate.id)\
                .filter(filter).all()
        return super(AuxService, self).prepare_external_values(**values)

    def after_create(self, new_object, **values):
        super(AuxService, self).after_create(new_object, **values)
        if values.get('template_id'):
            if isinstance(values['template_id'], (tuple, list, set)):
                ids = values.get('template_id')
            else:
                ids = [values.get('template_id')]
            for i in ids:
                get_or_create(self.db, AnnotationFormTemplateField,
                              form_field_id=new_object.id, form_template_id=i)
        return values

    def after_update(self, new_object, **values):
        super(AuxService, self).after_update(new_object, **values)
        if values.get('template_id'):
            if isinstance(values['template_id'], (tuple, list, set)):
                ids = values.get('template_id')
            else:
                ids = [values.get('template_id')]
            rl = self.db.query(AnnotationFormTemplateField)\
                .filter(AnnotationFormTemplateField.form_field_id == new_object.id)
            rl.filter(AnnotationFormTemplateField.form_template_id.notin_(ids)) \
                .delete(synchronize_session='fetch')
            for i in ids:
                get_or_create(self.db, AnnotationFormTemplateField,
                              form_field_id=new_object.id, form_template_id=i)
        return values

    def aux_filter(self, filter):
        clauses = []

        if filter.get('template_id'):
            _ids = self.db.query(AnnotationFormTemplateField.form_field_id) \
                .filter(filter_parse(AnnotationFormTemplateField, {'template_id': filter.get('template_id')}))
            clauses.append(self.orm.id.in_(_ids.subquery()))

        return clauses + super(AuxService, self).aux_filter(filter)