

def initialize_bibtex_forms():

	# Field types
	bibtex_fields = [
		{
			'name':
				'address',
			'description':
				"Publisher's address (usually just the city, but can be the full address for lesser-known publishers)",
		},
		{
			'name':
				'annote',
			'description':
				"An annotation for annotated bibliography styles (not typical)",
		},
		{
			'name':
				'author',
			'description':
				"The name(s) of the author(s) (in the case of more than one author, separated by and)",
		},
		{
			'name':
				'booktitle',
			'description':
				"The title of the book, if only part of it is being cited",
		},
		{
			'name':
				'Email',
			'description':
				"The email of the author(s)",
		},
		{
			'name':
				'chapter',
			'description':
				"The chapter number",
		},
		{
			'name':
				'crossref',
			'description':
				"The key of the cross-referenced entry",
		},
		{
			'name':
				'doi',
			'description':
				"digital object identifier",
		},
		{
			'name':
				'edition',
			'description':
				"The edition of a book, long form (such as \"First\" or \"Second\")",
		},
		{
			'name':
				'editor',
			'description':
				"The name(s) of the editor(s)",
		},
		{
			'name':
				'howpublished',
			'description':
				"How it was published, if the publishing method is nonstandard",
		},
		{
			'name':
				'institution',
			'description':
				"The institution that was involved in the publishing, but not necessarily the publisher",
		},
		{
			'name':
				'journal',
			'description':
				"The journal or magazine the work was published in",
		},
		{
			'name':
				'key',
			'description':
				"A hidden field used for specifying or overriding the alphabetical order of entries "
				"(when the \"author\" and \"editor\" fields are missing). Note that this is very different from the key "
				"(mentioned just after this list) that is used to cite or cross-reference the entry.",
		},
		{
			'name':
				'month',
			'description':
				"The month of publication (or, if unpublished, the month of creation)",
		},
		{
			'name':
				'note',
			'description':
				"Miscellaneous extra information",
		},
		{
			'name':
				'number',
			'description':
				"The \"(issue) number\" of a journal, magazine, or tech-report, if applicable. "
				"Note that this is not the \"article number\" assigned by some journals.",
		},
		{
			'name':
				'organization',
			'description':
				"The conference sponsor",
		},
		{
			'name':
				'pages',
			'description':
				"Page numbers, separated either by commas or double-hyphens.",
		},
		{
			'name':
				'publisher',
			'description':
				"The publisher's name",
		},
		{
			'name':
				'school',
			'description':
				"The school where the thesis was written",
		},
		{
			'name':
				'series',
			'description':
				"The series of books the book was published in "
				"(e.g. \"The Hardy Boys\" or \"Lecture Notes in Computer Science\")",
		},
		{
			'name':
				'title',
			'description':
				"The title of the work",
		},
		{
			'name':
				'type',
			'description':
				"The field overriding the default type of publication "
				"(e.g. \"Research Note\" for techreport, \"{PhD} dissertation\" for phdthesis, "
				"\"Section\" for inbook/incollection)",
		},
		{
			'name':
				'volume',
			'description':
				"The volume of a journal or multi-volume book",
		},
		{
			'name':
				'year',
			'description':
				"The year of publication (or, if unpublished, the year of creation)",
		},
	]
	from ..services.annotation_forms.fields import AuxService
	service = AuxService()
	fields = dict()
	try:
		for i in bibtex_fields:
			c = service.create(**i)[0]
			fields[c.name] = c
		service.db.commit()
	except Exception as e:
		print('Something went wrong when creating the BibTeX fields. They may already exist.')

	# Entry types / Templates

	bibtex_entities = [
		{
			'name':
				'article',
			'description':
				"An article from a journal or magazine.",
			'required_fields':
				['author', 'title', 'journal', 'year', 'volume'],
			'optional_fields':
				['number', 'pages', 'month', 'doi', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'book',
			'description':
				"A book with an explicit publisher.",
			'required_fields':
				['author', 'editor', 'title', 'publisher', 'year'],
			'optional_fields':
				['volume', 'number', 'series', 'address', 'edition', 'month', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'booklet',
			'description':
				"A work that is printed and bound, but without a named publisher or sponsoring institution.",
			'required_fields':
				['title'],
			'optional_fields':
				['author', 'howpublished', 'address', 'month', 'year', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'conference',
			'description':
				"The same as inproceedings, included for Scribe compatibility.",
			'required_fields':
				['author', 'title', 'booktitle', 'year'],
			'optional_fields':
				['editor', 'volume', 'number', 'series', 'pages', 'address', 'month', 'organization', 'publisher', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'inbook',
			'description':
				"A part of a book, usually untitled. May be a chapter (or section, etc.) and/or a range of pages.",
			'required_fields':
				['author', 'editor', 'title', 'chapter', 'pages', 'publisher', 'year'],
			'optional_fields':
				['volume', 'number', 'series', 'type', 'address', 'edition', 'month', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'incollection',
			'description':
				"A part of a book having its own title.",
			'required_fields':
				['author', 'title', 'booktitle', 'publisher', 'year'],
			'optional_fields':
				['editor', 'volume', 'number', 'series', 'type', 'chapter', 'pages', 'address', 'edition', 'month', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'inproceedings',
			'description':
				"An article in a conference proceedings.",
			'required_fields':
				['author', 'title', 'booktitle', 'year'],
			'optional_fields':
				['editor', 'volume', 'number', 'series', 'pages', 'address', 'month', 'organization', 'publisher', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'manual',
			'description':
				"Technical documentation.",
			'required_fields':
				['title'],
			'optional_fields':
				['author', 'organization', 'address', 'edition', 'month', 'year', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'mastersthesis',
			'description':
				"A master's thesis.",
			'required_fields':
				['author', 'title', 'school', 'year'],
			'optional_fields':
				['type', 'address', 'month', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'misc',
			'description':
				"For use when nothing else fits.",
			'required_fields':
				[],
			'optional_fields':
				['author', 'title', 'howpublished', 'month', 'year', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'phdthesis',
			'description':
				"A Ph.D. thesis.",
			'required_fields':
				['author', 'title', 'school', 'year'],
			'optional_fields':
				['type', 'address', 'month', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'proceedings',
			'description':
				"The proceedings of a conference.",
			'required_fields':
				['title', 'year'],
			'optional_fields':
				['editor', 'volume', 'number', 'series', 'address', 'month', 'publisher', 'organization', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'techreport',
			'description':
				"A report published by a school or other institution, usually numbered within a series.",
			'required_fields':
				['author', 'title', 'institution', 'year'],
			'optional_fields':
				['type', 'number', 'address', 'month', 'note', 'key'],
			'type':
				'template_bibtex',
		},
		{
			'name':
				'unpublished',
			'description':
				"A document having an author and title, but not formally published.",
			'required_fields':
				['author', 'title', 'note'],
			'optional_fields':
				['month', 'year', 'key'],
			'type':
				'template_bibtex',
		},
	]
	from ..services.annotation_forms.templates import AuxService
	service = AuxService()
	try:
		for i in bibtex_entities:
			i['required_fields'] = [fields.get(name, name) for name in i.get('required_fields', [])]
			i['optional_fields'] = [fields.get(name, name) for name in i.get('optional_fields', [])]
			service.create(**i)
		service.db.commit()
	except Exception as e:
		print('Something went wrong when creating the BibTeX templates. They may already exist.')
