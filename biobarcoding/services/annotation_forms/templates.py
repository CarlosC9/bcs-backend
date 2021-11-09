from . import FormItemAuxService
from ..main import get_orm
from ...rest import filter_parse
from ...db_models import DBSession
from ...db_models.sysadmin import AnnotationFormTemplateField


##
# TEMPLATE TOOLS
##

class AuxService(FormItemAuxService):

    def __init__(self):
        self.orm = get_orm('templates')

    def aux_filter(self, filter):
        clauses = []

        if filter.get('field_id'):
            _ids = DBSession.query(AnnotationFormTemplateField.form_template_id) \
                .filter(filter_parse(AnnotationFormTemplateField, {'field_id': filter.get('field_id')}))
            clauses.append(self.orm.id.in_(_ids))   # .all()

        return clauses + super(AuxService, self).aux_filter()