from ..main import BasicService, get_orm


##
# ANNOTATION FORM SERVICE
##

class FormRelationshipService(BasicService):

	def __init__(self):
		super(FormRelationshipService, self).__init__()
		self.orm = get_orm('form_relationship')

	def prepare_values(self, template=None, field=None, **values):

		form_template = values.get('form_template')
		if isinstance(form_template, str):
			template = template or form_template
			values.pop('form_template')
		if template and not values.get('form_template') and not values.get('form_template_id'):
			from ...db_models.sa_annotations import AnnotationFormTemplate
			values['form_template'] = self.db.query(AnnotationFormTemplate) \
				.filter(AnnotationFormTemplate.name == template).one()

		form_field = values.get('form_field')
		if isinstance(form_field, str):
			field = field or form_field
			values.pop('form_field')
		if field and not values.get('form_field') and not values.get('form_field_id'):
			from ...db_models.sa_annotations import AnnotationFormField
			values['form_field'] = self.db.query(AnnotationFormField) \
				.filter(AnnotationFormField.name == field).one()

		return super(FormRelationshipService, self).prepare_values(**values)


##
# ANNOTATION FORM SERVICE
##

class RelationshipService(BasicService):

	def __init__(self):
		super(RelationshipService, self).__init__()
		self.orm = get_orm('relationship')

	def prepare_values(self, rl=None, **values):

		subject_uuid = values.get('subject_uuid')
		if subject_uuid and not values.get('subject_id') and not values.get('subject'):
			from ...db_models.bioinformatics import FunctionalObject
			values['subject'] = self.db.query(FunctionalObject) \
				.filter(FunctionalObject.uuid == subject_uuid).one()

		object_uuid = values.get('object_uuid')
		if object_uuid and not values.get('object_id') and not values.get('object'):
			from ...db_models.bioinformatics import FunctionalObject
			values['object'] = self.db.query(FunctionalObject) \
				.filter(FunctionalObject.uuid == object_uuid).one()

		# form_rl = values.get('form_rl')
		# if isinstance(form_rl, str):
		# 	rl = rl or form_rl
		# 	values.pop('form_rl')
		# if rl and not values.get('form_rl') and not values.get('form_rl_id'):
		# 	from ...db_models.sa_annotations import AnnotationFormTemplateField
		# 	values['form_rl'] = self.db.query(AnnotationFormTemplateField) \
		# 		.filter(AnnotationFormTemplateField.name == rl).one()

		return super(RelationshipService, self).prepare_values(**values)
