import json
from urllib.error import URLError

from .. import REQUEST_URL
from ...system import SA_TASK_SESSION
from ....services import log_exception


def _create_request(service, **kwargs):
	try:
		url = f"{REQUEST_URL}/annotation_form_{service}/"
		print('POST ' + url)
		return SA_TASK_SESSION.post(url, json=kwargs, headers={'Content-Type': 'application/json'})
	except Exception as e:
		print(f'Something went wrong when creating the {kwargs.get("name")}. It may already exist.')
		log_exception(e)
		return None


object_type_id = []


##
# BibTex
##

def initialize_bibtex_forms():

	# TODO: create cv and cvterms

	print(' > Creating BibTex fields')
	bibtex_fields = [
		{
			'name':
				'address',
			'cvterm':
				'address',
			'description':
				"Publisher's address (usually just the city, but can be the full address for lesser-known publishers)",
		},
		{
			'name':
				'annote',
			'cvterm':
				'annote',
			'description':
				"An annotation for annotated bibliography styles (not typical)",
		},
		{
			'name':
				'author',
			'cvterm':
				'author',
			'description':
				"The name(s) of the author(s) (in the case of more than one author, separated by and)",
		},
		{
			'name':
				'booktitle',
			'cvterm':
				'booktitle',
			'description':
				"The title of the book, if only part of it is being cited",
		},
		{
			'name':
				'Email',
			'cvterm':
				'Email',
			'description':
				"The email of the author(s)",
		},
		{
			'name':
				'chapter',
			'cvterm':
				'chapter',
			'description':
				"The chapter number",
		},
		{
			'name':
				'crossref',
			'cvterm':
				'crossref',
			'description':
				"The key of the cross-referenced entry",
		},
		{
			'name':
				'doi',
			'cvterm':
				'doi',
			'description':
				"digital object identifier",
		},
		{
			'name':
				'edition',
			'cvterm':
				'edition',
			'description':
				"The edition of a book, long form (such as \"First\" or \"Second\")",
		},
		{
			'name':
				'editor',
			'cvterm':
				'editor',
			'description':
				"The name(s) of the editor(s)",
		},
		{
			'name':
				'howpublished',
			'cvterm':
				'howpublished',
			'description':
				"How it was published, if the publishing method is nonstandard",
		},
		{
			'name':
				'institution',
			'cvterm':
				'institution',
			'description':
				"The institution that was involved in the publishing, but not necessarily the publisher",
		},
		{
			'name':
				'journal',
			'cvterm':
				'journal',
			'description':
				"The journal or magazine the work was published in",
		},
		{
			'name':
				'key',
			'cvterm':
				'key',
			'description':
				"A hidden field used for specifying or overriding the alphabetical order of entries "
				"(when the \"author\" and \"editor\" fields are missing). Note that this is very different from the key "
				"(mentioned just after this list) that is used to cite or cross-reference the entry.",
		},
		{
			'name':
				'month',
			'cvterm':
				'month',
			'description':
				"The month of publication (or, if unpublished, the month of creation)",
		},
		{
			'name':
				'note',
			'cvterm':
				'note',
			'description':
				"Miscellaneous extra information",
		},
		{
			'name':
				'number',
			'cvterm':
				'number',
			'description':
				"The \"(issue) number\" of a journal, magazine, or tech-report, if applicable. "
				"Note that this is not the \"article number\" assigned by some journals.",
		},
		{
			'name':
				'organization',
			'cvterm':
				'organization',
			'description':
				"The conference sponsor",
		},
		{
			'name':
				'pages',
			'cvterm':
				'pages',
			'description':
				"Page numbers, separated either by commas or double-hyphens.",
		},
		{
			'name':
				'publisher',
			'cvterm':
				'publisher',
			'description':
				"The publisher's name",
		},
		{
			'name':
				'school',
			'cvterm':
				'school',
			'description':
				"The school where the thesis was written",
		},
		{
			'name':
				'series',
			'cvterm':
				'series',
			'description':
				"The series of books the book was published in "
				"(e.g. \"The Hardy Boys\" or \"Lecture Notes in Computer Science\")",
		},
		{
			'name':
				'title',
			'cvterm':
				'title',
			'description':
				"The title of the work",
		},
		{
			'name':
				'type',
			'cvterm':
				'type',
			'description':
				"The field overriding the default type of publication "
				"(e.g. \"Research Note\" for techreport, \"{PhD} dissertation\" for phdthesis, "
				"\"Section\" for inbook/incollection)",
		},
		{
			'name':
				'volume',
			'cvterm':
				'volume',
			'description':
				"The volume of a journal or multi-volume book",
		},
		{
			'name':
				'year',
			'cvterm':
				'year',
			'description':
				"The year of publication (or, if unpublished, the year of creation)",
		},
	]
	for i in bibtex_fields:
		_create_request('fields', **i, cv='bibtex', standard='BibTex', object_type_id=object_type_id)

	print(' > Creating BibTex templates')
	bibtex_entities = [
		{
			'name':
				'article',
			'cvterm':
				'article',
			'description':
				"An article from a journal or magazine.",
			'required_field':
				['author', 'title', 'journal', 'year', 'volume'],
			'field':
				['number', 'pages', 'month', 'doi', 'note', 'key'],
		},
		{
			'name':
				'book',
			'cvterm':
				'book',
			'description':
				"A book with an explicit publisher.",
			'required_field':
				['author', 'editor', 'title', 'publisher', 'year'],
			'field':
				['volume', 'number', 'series', 'address', 'edition', 'month', 'note', 'key'],
		},
		{
			'name':
				'booklet',
			'cvterm':
				'booklet',
			'description':
				"A work that is printed and bound, but without a named publisher or sponsoring institution.",
			'required_field':
				['title'],
			'field':
				['author', 'howpublished', 'address', 'month', 'year', 'note', 'key'],
		},
		{
			'name':
				'conference',
			'cvterm':
				'conference',
			'description':
				"The same as inproceedings, included for Scribe compatibility.",
			'required_field':
				['author', 'title', 'booktitle', 'year'],
			'field':
				['editor', 'volume', 'number', 'series', 'pages', 'address', 'month', 'organization', 'publisher', 'note', 'key'],
		},
		{
			'name':
				'inbook',
			'cvterm':
				'inbook',
			'description':
				"A part of a book, usually untitled. May be a chapter (or section, etc.) and/or a range of pages.",
			'required_field':
				['author', 'editor', 'title', 'chapter', 'pages', 'publisher', 'year'],
			'field':
				['volume', 'number', 'series', 'type', 'address', 'edition', 'month', 'note', 'key'],
		},
		{
			'name':
				'incollection',
			'cvterm':
				'incollection',
			'description':
				"A part of a book having its own title.",
			'required_field':
				['author', 'title', 'booktitle', 'publisher', 'year'],
			'field':
				['editor', 'volume', 'number', 'series', 'type', 'chapter', 'pages', 'address', 'edition', 'month', 'note', 'key'],
		},
		{
			'name':
				'inproceedings',
			'cvterm':
				'inproceedings',
			'description':
				"An article in a conference proceedings.",
			'required_field':
				['author', 'title', 'booktitle', 'year'],
			'field':
				['editor', 'volume', 'number', 'series', 'pages', 'address', 'month', 'organization', 'publisher', 'note', 'key'],
		},
		{
			'name':
				'manual',
			'cvterm':
				'manual',
			'description':
				"Technical documentation.",
			'required_field':
				['title'],
			'field':
				['author', 'organization', 'address', 'edition', 'month', 'year', 'note', 'key'],
		},
		{
			'name':
				'mastersthesis',
			'cvterm':
				'mastersthesis',
			'description':
				"A master's thesis.",
			'required_field':
				['author', 'title', 'school', 'year'],
			'field':
				['type', 'address', 'month', 'note', 'key'],
		},
		{
			'name':
				'misc',
			'cvterm':
				'misc',
			'description':
				"For use when nothing else fits.",
			'required_field':
				[],
			'field':
				['author', 'title', 'howpublished', 'month', 'year', 'note', 'key'],
		},
		{
			'name':
				'phdthesis',
			'cvterm':
				'phdthesis',
			'description':
				"A Ph.D. thesis.",
			'required_field':
				['author', 'title', 'school', 'year'],
			'field':
				['type', 'address', 'month', 'note', 'key'],
		},
		{
			'name':
				'proceedings',
			'cvterm':
				'proceedings',
			'description':
				"The proceedings of a conference.",
			'required_field':
				['title', 'year'],
			'field':
				['editor', 'volume', 'number', 'series', 'address', 'month', 'publisher', 'organization', 'note', 'key'],
		},
		{
			'name':
				'techreport',
			'cvterm':
				'techreport',
			'description':
				"A report published by a school or other institution, usually numbered within a series.",
			'required_field':
				['author', 'title', 'institution', 'year'],
			'field':
				['type', 'number', 'address', 'month', 'note', 'key'],
		},
		{
			'name':
				'unpublished',
			'cvterm':
				'unpublished',
			'description':
				"A document having an author and title, but not formally published.",
			'required_field':
				['author', 'title', 'note'],
			'field':
				['month', 'year', 'key'],
		},
	]
	for i in bibtex_entities:
		_create_request('templates', **i, cv='bibtex', standard='BibTex', object_type_id=object_type_id)


##
# Simple Darwin Core
##

def __dwc_term2annotation(**kwargs):
	# TODO: 'term_modified', 'term_deprecated', 'replaces_term', 'replaces1_term'
	ann = {
		'cv': kwargs.get('database'),		# (p.e. 'terms')
		'db': kwargs.get('vann_preferredNamespacePrefix'),		# (p.e. 'dwc')
		'standard': 'Darwin Core',		# kwargs.get('standard'),		# (p.e. 'http://www.tdwg.org/standards/450')
		'cvterm': kwargs.get('term_localName'),		# (p.e. 'acceptedNameUsage')
		'url': kwargs.get('vann_preferredNamespaceUri', '') + kwargs.get('term_localName', ''),
			# (p.e. 'http://rs.tdwg.org/dwc/terms/acceptedScientificName')
		'name': kwargs.get('label'),		# (p.e. 'Accepted Name Usage')
		'rdf_type': kwargs.get('rdf_type'),		# (p.e. 'http://www.w3.org/1999/02/22-rdf-syntax-ns#Property')
		'description': kwargs.get('rdfs_comment'),		# (p.e. )
		'definition': kwargs.get('dcterms_description'),		# (p.e. )
	}
	if kwargs.get('tdwgutility_organizedInClass'):
		ann['template'] = kwargs.get('tdwgutility_organizedInClass').rsplit('/', 1)[-1]
	ann['object_type_id'] = object_type_id
	return ann


def initialize_dwc_forms():
	# EDITED:
	# Script to build Markdown pages that provide term metadata for complex vocabularies
	# Steve Baskauf 2020-08-12 CC0
	# This script merges static Markdown header and footer documents with term information tables (in Markdown)
	# generated from data in the rs.tdwg.org repo from the TDWG Github site

	import pandas as pd

	# This is the base URL for raw files from the branch of the repo that has been pushed to GitHub
	githubBaseUri = 'https://raw.githubusercontent.com/tdwg/rs.tdwg.org/master/'

	# This is a Python list of the database names of the term lists to be included in the document.
	terms_dbb = ['terms', 'iri', 'dc-for-dwc', 'dcterms-for-dwc']
				# + ['audubon', 'curatorial', 'dwc-obsolete', 'dwcore', 'dwctype', 'geospatial', 'utility']

	print('\n > Retrieving namespaces from GitHub')
	n, namespaces = set(), []

	try:
		frame = pd.read_csv(githubBaseUri + 'term-lists/term-lists.csv', na_filter=False)
	except Exception as e:
		print('\n > ERROR: The DwC data could not be retrieved')
		return None

	for db in terms_dbb:
		for index, row in frame.iterrows():
			if row.get('database') == db:
				namespaces.append(row)
				break
			elif row.get('database') not in terms_dbb:
				n.add(row.get('database'))
	print(f'NOT INCLUDED: {len(n)}\n', sorted(n))
	print(f'INCLUDED: {len(namespaces)}\n', [i.get('database') for i in namespaces])

	print('\n > Retrieving metadata of terms from all selected namespaces from GitHub')

	# Create column list
	column_list = [
		'vann_preferredNamespacePrefix', 'vann_preferredNamespaceUri', 'term_localName', 'label', 'rdfs_comment',
		'dcterms_description', 'examples', 'term_modified', 'term_deprecated', 'rdf_type', 'replaces_term',
		'replaces1_term', 'controlled_value_string', 'tdwgutility_organizedInClass']
	# + ['version_iri']

	# Create list of lists metadata table
	term_lists = []
	for namespace in namespaces:
		# TODO: create cv / db ?
		try:
			# retrieve current term metadata for term list
			data_url = githubBaseUri + namespace.get('database', '') + '/' + namespace.get('database', '') + '.csv'
			frame = pd.read_csv(data_url, na_filter=False)
		except URLError as e:
			print(namespace.get('database') + ' could not be found')
			continue
		for index, row in frame.iterrows():
			term_list = [
				namespace['vann_preferredNamespacePrefix'], namespace['vann_preferredNamespaceUri'],
				row['term_localName'], row['label'], row['rdfs_comment'],
				row['dcterms_description'], row['examples'], row['term_modified'], row['term_deprecated'],
				row['rdf_type'], row['replaces_term'], row['replaces1_term'],
				row.get('controlled_value_string'), row['tdwgutility_organizedInClass']]
			# TODO: create cvterm ?
			term_lists.append(term_list)
	print(f'INCLUDED TERMS: {len(term_lists)}')

	# Turn list of lists into dataframe
	terms_df = pd.DataFrame(term_lists, columns=column_list)
	terms_df = terms_df[terms_df.term_deprecated != 'true']		# TODO: store and manage deprecated terms ?

	# This makes sort case insensitive
	terms_sorted_by_localname = terms_df.iloc[terms_df.term_localName.str.lower().argsort()]

	# # To organize in categories, the display_order list must contain the IRIs that are values of tdwgutility_organizedInClass
	# display_order = ['', 'http://purl.org/dc/elements/1.1/', 'http://purl.org/dc/terms/', 'http://rs.tdwg.org/dwc/terms/attributes/UseWithIRI']
	# display_label = ['Record level', 'Dublin Core legacy namespace', 'Dublin Core terms namespace', 'IRI-value terms']
	category_order = list(set(terms_df.tdwgutility_organizedInClass))
	_ = [i.rsplit('/', 1) for i in category_order]
	category_label = [i[-1] for i in _]
	_ = [(i[0]+'/', i[-1]) for i in _]
	groups_df = terms_df[terms_df[['vann_preferredNamespaceUri', 'term_localName']].apply(tuple, axis=1).isin(_)]
	groups = set(category_label)

	print('\n > Creating Darwin Core templates')

	print('**Classes**')
	classes = set(terms_df[terms_df.rdf_type.str.endswith('rdf-schema#Class')]['term_localName'])
	templates = classes.intersection(groups)
	# TODO: what to do?
	#  not_templates = classes.difference(groups)
	#  free_fields = groups.difference(classes)
	for row_index, row in groups_df.iterrows():
		print('creating template ' + row.get('label'))
		uri = row['vann_preferredNamespaceUri'] + row['term_localName']
		curie = row['vann_preferredNamespacePrefix'] + ":" + row['term_localName']
		print('[' + curie + '](#' + curie.replace(':', '_') + ')\n' + uri)
		_create_request('templates', **__dwc_term2annotation(**row))

	print('\n > Creating Darwin Core fields')

	for i in range(0, len(category_order)):
		print('\n**' + category_label[i] + '** (' + category_order[i] + ')')
		filtered_table = terms_sorted_by_localname[terms_sorted_by_localname['tdwgutility_organizedInClass'] == category_order[i]]
		filtered_table.reset_index(drop=True, inplace=True)

		for row_index, row in filtered_table.iterrows():
			# print('creating field', row.get('label'))
			# uri = row['vann_preferredNamespaceUri'] + row['term_localName']
			# curie = row['vann_preferredNamespacePrefix'] + ":" + row['term_localName']
			# print('[' + curie + '](#' + curie.replace(':', '_') + ')\n' + uri)
			values = __dwc_term2annotation(**row)
			if category_label[i] not in templates:
				values.pop('tdwgutility_organizedInClass', '')
			_create_request('fields', **values)


def run():
	url = f"{REQUEST_URL}/object_types/"
	global object_type_id
	print('GET ' + url)
	object_type_id = [i['id'] for i in json.loads(SA_TASK_SESSION.get(url).text)['content']]
	initialize_bibtex_forms()
	initialize_dwc_forms()
	return 'DONE'
