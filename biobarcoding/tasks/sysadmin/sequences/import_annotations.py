import json

from biobarcoding.services import log_exception
from biobarcoding.tasks.sysadmin import create_request, read_request, update_request


def run(ann_entries: dict):

	def load_response(r):
		try:
			return json.loads(r.text)
		except Exception as e:
			return {}

	def get_response_count(r):
		return load_response(r).get('count')

	def get_response_content(r):
		return load_response(r).get('content')

	def get_response_id(r):
		try:
			c = get_response_content(r)
			if isinstance(c, (tuple, list, set)):
				c = c[0] if len(c) == 1 else {}
			return c.get('id')
		except Exception as e:
			print('missing id for ', c)
			return None

	def get_ann_type(arg):
		return'template' if arg.get('template') else 'field' if arg.get('field') else None

	def build_template(template_id, standard, value):
		for field_name in value.copy():
			try:
				_ = read_request('annotation_form_fields/',
									name=field_name, template_id=template_id, standard=standard)
				if get_response_count(_) != 1:
					raise Exception()
			except Exception as e:
				_ = create_request('annotation_form_fields/',
									name=field_name, template_id=template_id, standard=standard)
			if get_response_id(_):
				value[get_response_id(_)] = value.pop(field_name)
		print(list(value.keys()))
		return value

	c_form, c_ann = 0, 0
	for _hash, seqs in ann_entries.items():

		ann = json.loads(_hash)
		form_item_type = get_ann_type(ann)
		if not form_item_type:
			continue

		form_item = ann.get(form_item_type)
		form_item = form_item if isinstance(form_item, dict) else json.loads(form_item)
		form_value = ann.get('value')
		try:
			form_json_value = json.loads(form_value)
		except Exception as e:
			form_json_value = form_value

		# get_or_create form_item
		try:
			response = read_request('annotation_form_%ss/' % form_item_type, filter=form_item)
			if get_response_count(response) != 1:
				raise Exception()
		except Exception as e:
			print('"form_%s" not found, creating it' % form_item_type)
			response = create_request('annotation_form_%ss/' % form_item_type, **form_item)
			c_form += 1

		if form_item_type == 'template':
			form_template = get_response_content(response)[0]
			build_template(form_template.get('id'), form_template.get('standard'), form_json_value)

		# get_or_create template or field
		form_item_id = get_response_id(response)
		if not form_item_id:
			continue
		values = {'form_%s_id' % form_item_type: form_item_id, 'type': form_item_type}
		try:
			response = read_request('annotations', filter={**values, 'value': form_value})
			if get_response_count(response) != 1:
				raise Exception()
		except Exception as e:
			print('"%s_item" not found, creating it' % form_item_type)
			create_request('annotations', **values, object_uuid=seqs, value=form_json_value)
			c_ann += 1
			continue

		# put annotation links
		ann_id = get_response_id(response)
		try:
			update_request('annotations/%s' % ann_id, new_object_uuid=seqs)
		except Exception as e:
			log_exception(e)
			return 'EXCEPTION'

	print('ANNOTATION_FORMS:', c_form)
	print('ANNOTATION_ITEMS:', c_ann)
	print('ANNOTATION_LINKS:', len([id for seq_ids in ann_entries.values() for id in seq_ids]))

	return 'DONE'
