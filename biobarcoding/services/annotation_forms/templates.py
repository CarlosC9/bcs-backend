from . import FormItemAuxService
from ..main import get_orm
from ...rest import filter_parse
from ...services import get_or_create
from ...db_models import DBSession
from ...db_models.sa_annotations import AnnotationFormField, AnnotationFormTemplateField


##
# TEMPLATE TOOLS
##

class AuxService(FormItemAuxService):

    def __init__(self):
        super(AuxService, self).__init__()
        self.orm = get_orm('templates')

    def after_create(self, new_object, field=[], **kwargs):
        if field and not kwargs.get('field_id'):
            if isinstance(field, (tuple, list, set)):
                filter = AnnotationFormField.name.in_(field)
            else:
                filter = AnnotationFormField.name == field
            kwargs['field_id'] = DBSession.query(AnnotationFormField.id)\
                .filter(filter).all()
        if kwargs.get('field_id'):
            if isinstance(kwargs['field_id'], (tuple, list, set)):
                ids = kwargs.get('field_id')
            else:
                ids = [kwargs.get('field_id')]
            for i in ids:
                get_or_create(DBSession, AnnotationFormTemplateField,
                              form_template_id=new_object.id, form_field_id=i)
        return super(AuxService, self).after_create(new_object, field=field, **kwargs)

    def aux_filter(self, filter):
        clauses = []

        if filter.get('field_id'):
            _ids = DBSession.query(AnnotationFormTemplateField.form_template_id) \
                .filter(filter_parse(AnnotationFormTemplateField, {'field_id': filter.get('field_id')}))
            clauses.append(self.orm.id.in_(_ids.subquery()))

        return clauses + super(AuxService, self).aux_filter(filter)