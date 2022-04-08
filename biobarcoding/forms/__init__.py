

def initialize_bibtex_forms():

	from ..db_models import DBSession, ObjectType
	object_type_id = [i.id for i in DBSession.query(ObjectType).all()]

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
	from ..services.annotation_forms.fields import Service
	service = Service()
	fields = dict()
	for i in bibtex_fields:
		try:
			c = service.create(**i, object_type_id=object_type_id)[0]
			fields[c.name] = c
			service.db.commit()
		except:
			service.db.rollback()
			print(f'Something went wrong when creating the {i.get("name")} BibTeX field. It may already exist.')

	# Entry types / Templates

	bibtex_entities = [
		{
			'name':
				'article',
			'description':
				"An article from a journal or magazine.",
			'required_field':
				['author', 'title', 'journal', 'year', 'volume'],
			'field':
				['number', 'pages', 'month', 'doi', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'book',
			'description':
				"A book with an explicit publisher.",
			'required_field':
				['author', 'editor', 'title', 'publisher', 'year'],
			'field':
				['volume', 'number', 'series', 'address', 'edition', 'month', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'booklet',
			'description':
				"A work that is printed and bound, but without a named publisher or sponsoring institution.",
			'required_field':
				['title'],
			'field':
				['author', 'howpublished', 'address', 'month', 'year', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'conference',
			'description':
				"The same as inproceedings, included for Scribe compatibility.",
			'required_field':
				['author', 'title', 'booktitle', 'year'],
			'field':
				['editor', 'volume', 'number', 'series', 'pages', 'address', 'month', 'organization', 'publisher', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'inbook',
			'description':
				"A part of a book, usually untitled. May be a chapter (or section, etc.) and/or a range of pages.",
			'required_field':
				['author', 'editor', 'title', 'chapter', 'pages', 'publisher', 'year'],
			'field':
				['volume', 'number', 'series', 'type', 'address', 'edition', 'month', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'incollection',
			'description':
				"A part of a book having its own title.",
			'required_field':
				['author', 'title', 'booktitle', 'publisher', 'year'],
			'field':
				['editor', 'volume', 'number', 'series', 'type', 'chapter', 'pages', 'address', 'edition', 'month', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'inproceedings',
			'description':
				"An article in a conference proceedings.",
			'required_field':
				['author', 'title', 'booktitle', 'year'],
			'field':
				['editor', 'volume', 'number', 'series', 'pages', 'address', 'month', 'organization', 'publisher', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'manual',
			'description':
				"Technical documentation.",
			'required_field':
				['title'],
			'field':
				['author', 'organization', 'address', 'edition', 'month', 'year', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'mastersthesis',
			'description':
				"A master's thesis.",
			'required_field':
				['author', 'title', 'school', 'year'],
			'field':
				['type', 'address', 'month', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'misc',
			'description':
				"For use when nothing else fits.",
			'required_field':
				[],
			'field':
				['author', 'title', 'howpublished', 'month', 'year', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'phdthesis',
			'description':
				"A Ph.D. thesis.",
			'required_field':
				['author', 'title', 'school', 'year'],
			'field':
				['type', 'address', 'month', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'proceedings',
			'description':
				"The proceedings of a conference.",
			'required_field':
				['title', 'year'],
			'field':
				['editor', 'volume', 'number', 'series', 'address', 'month', 'publisher', 'organization', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'techreport',
			'description':
				"A report published by a school or other institution, usually numbered within a series.",
			'required_field':
				['author', 'title', 'institution', 'year'],
			'field':
				['type', 'number', 'address', 'month', 'note', 'key'],
			'standard':
				'bibtex',
		},
		{
			'name':
				'unpublished',
			'description':
				"A document having an author and title, but not formally published.",
			'required_field':
				['author', 'title', 'note'],
			'field':
				['month', 'year', 'key'],
			'standard':
				'bibtex',
		},
	]
	from ..services.annotation_forms.templates import Service
	service = Service()
	for i in bibtex_entities:
		try:
			service.create(**i, object_type_id=object_type_id)
			service.db.commit()
		except:
			service.db.rollback()
			print(f'Something went wrong when creating the {i.get("name")} BibTeX template. It may already exist.')
