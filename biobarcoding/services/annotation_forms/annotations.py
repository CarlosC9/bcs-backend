from ..main import SimpleAuxService, get_orm


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
                from ...db_models.sysadmin import AnnotationFormTemplate
                values['form_template'] = self.db.query(AnnotationFormTemplate)\
                    .filter(AnnotationFormTemplate.name == template).one()
        if field:
            if not template and not values.get("type"):
                values["type"] = 'field'
            if not values.get('form_field') and not values.get('form_field_id'):
                from ...db_models.sysadmin import AnnotationFormField
                values['form_field'] = self.db.query(AnnotationFormField)\
                    .filter(AnnotationFormField.name == field).one()

        if values.get("type") in ('template', 'field', 'text'):
            self.orm = get_orm(f'annotation_{values.get("type")}')

        return super(AuxService, self).prepare_values(**values)
