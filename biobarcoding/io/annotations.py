import json

from ..db_models import DBSession
from ..db_models.core import data_object_type_id
from ..db_models.sa_annotations import AnnotationFormTemplate, AnnotationFormField, AnnotationFormItemObjectType, \
	AnnotationFormTemplateField, AnnotationTemplate, AnnotationField
from ..services import get_or_create

BATCH_SIZE = 500
ANN_ENTRIES = {}


def clear_entries():
	ANN_ENTRIES.clear()


def get_or_create_ann_form(name: str, standard: str = '', fields: list = (), object_types: list = ()):
	form_template = get_or_create(DBSession, AnnotationFormTemplate, no_flush=True, name=name)
	if not form_template.id:
		DBSession.add(form_template)
		form_template.standard = form_template.standard or standard
		for obj_type in object_types:
			get_or_create(DBSession, AnnotationFormItemObjectType, no_flush=True,
						form_item=form_template, object_type_id=data_object_type_id[obj_type])
		for field in fields:
			form_field = get_or_create(DBSession, AnnotationFormField, no_flush=True, name=field)
			DBSession.add(form_field)
			form_field.standard = form_field.standard or standard
			get_or_create(DBSession, AnnotationFormItemObjectType, no_flush=True,
						  form_item=form_field, object_type_id=data_object_type_id['sequence'])
			get_or_create(DBSession, AnnotationFormTemplateField, no_flush=True,
						  form_template=form_template, form_field=form_field)
	return form_template


def ann_value_dump(arg):
	try:
		return json.dumps(arg, sort_keys=True)
	except:
		if isinstance(arg, (tuple, list, set)):
			return [ann_value_dump(_) for _ in arg]
		try:
			return ann_value_dump(arg.__dict__)
		except:
			try:
				return ann_value_dump(dict(((k, ann_value_dump(v)) for k, v in arg.items())))
			except:
				try:
					return ann_value_dump(dict(((k, ann_value_dump(v)) for k, v in arg.__dict__.items())))
				except:
					return str(arg)


def get_or_create_ann_template(value: any, standard: str = '', **kwargs):

	def encode_template_value(v: dict, s: str = ''):
		f_value = v.copy()
		for k, v in v.items():
			form_field = get_or_create(DBSession, AnnotationFormField, no_flush=True, standard=s, name=k)
			f_value[form_field.id] = f_value.get(form_field.id, f_value.pop(form_field.name, None))
		return f_value

	# get_or_create for value(JSONB)
	form_value = encode_template_value(value if isinstance(value, dict) else json.loads(ann_value_dump(value)), standard)
	instance = DBSession.query(AnnotationTemplate).filter_by(value=ann_value_dump(form_value), **kwargs).first()
	if not instance:
		instance = AnnotationTemplate(value=json.loads(ann_value_dump(form_value)), **kwargs)
	return instance


def import_file():
	pass
