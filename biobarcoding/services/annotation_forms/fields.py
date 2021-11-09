from . import FormItemAuxService
from ..main import get_orm
from ...rest import filter_parse
from ...services import get_or_create
from ...db_models import DBSession
from ...db_models.sysadmin import AnnotationFormTemplate, AnnotationFormTemplateField


##
# FIELD TOOLS
##

class AuxService(FormItemAuxService):

    def __init__(self):
        self.orm = get_orm('fields')

    def after_create(self, new_object, template=[], **kwargs):
        if template and not kwargs.get('template_id'):
            if hasattr(template, '__iter__'):
                filter = AnnotationFormTemplate.name.in_(template)
            else:
                filter = AnnotationFormTemplate.name == template
            kwargs['template_id'] = DBSession.query(AnnotationFormTemplate.id)\
                .filter(filter).all()
        if kwargs.get('template_id'):
            if hasattr(kwargs['template_id'], '__iter__'):
                ids = kwargs.get('template_id')
            else:
                ids = [kwargs.get('template_id')]
            for i in ids:
                get_or_create(DBSession, AnnotationFormTemplateField,
                              form_field_id=new_object.id, form_template_id=i)
        return super(AuxService, self).after_create()

    def aux_filter(self, filter):
        clauses = []

        if filter.get('template_id'):
            _ids = DBSession.query(AnnotationFormTemplateField.form_field_id) \
                .filter(filter_parse(AnnotationFormTemplateField, {'template_id': filter.get('template_id')}))
            clauses.append(self.orm.id.in_(_ids))   # .all()

        return clauses + super(AuxService, self).aux_filter()