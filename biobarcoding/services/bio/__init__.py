from sqlalchemy import inspect

from .. import get_bioformat
from ..main import BasicService
from ...db_models import DBSessionChado
from ...rest import filter_parse


##
# CHADO TOOLS
##

class BioService(BasicService):

	def __init__(self):
		super(BioService, self).__init__()
		self.db = DBSessionChado

	def prepare_values(self, cv=None, cvterm=None, db=None, dbxref=None, **values):

		from ...db_models import DBSessionChado
		from ...db_models.chado import Cv, Cvterm, Db, Dbxref

		if cv and cvterm and not values.get('cvterm_id'):
			values['cvterm_id'] = DBSessionChado.query(Cvterm).join(Cv) \
				.filter(Cv.name == cv, Cvterm.name == cvterm).one().cvterm_id

		if db and dbxref and not values.get('dbxref_id'):
			values['dbxref_id'] = DBSessionChado.query(Dbxref).join(Db) \
				.filter(Db.name == db, Dbxref.accession == dbxref).one().dbxref_id

		if values.get('cvterm_id') and not values.get('dbxref_id'):
			values['dbxref_id'] = DBSessionChado.query(Cvterm) \
				.filter(Cvterm.cvterm_id == values.get('cvterm_id')) \
				.first().dbxref_id

		if values.get('dbxref_id') and not values.get('cvterm_id'):
			values['cvterm_id'] = DBSessionChado.query(Cvterm) \
				.filter(Cvterm.dbxref_id == values.get('dbxref_id')) \
				.first().cvterm_id

		return super(BioService, self).prepare_values(**values)

	def aux_filter(self, _filter: dict) -> list:
		clauses = []

		_v = _filter.get('transaction_id') or _filter.get('version') or _filter.get('version_id')
		if _v:
			from sqlalchemy_continuum import version_class
			Version = version_class(self.orm)
			_ids = self.db.query(inspect(Version).primary_key[0]).filter(filter_parse(Version, {'transaction_id': _v}))
			clauses.append(inspect(self.orm).primary_key[0].in_(_ids))

		if _filter.get('issued_at'):
			from sqlalchemy_continuum import transaction_class
			Transaction = transaction_class(self.orm)
			_ids = self.db.query(Transaction.id).filter(filter_parse(Transaction, {'issued_at': _filter.get('issued_at')}))
			from sqlalchemy_continuum import version_class
			Version = version_class(self.orm)
			_ids = self.db.query(inspect(Version).primary_key[0]).filter(Version.transaction_id.in_(_ids))
			clauses.append(inspect(self.orm).primary_key[0].in_(_ids))

		return clauses + super(BioService, self).aux_filter(_filter)

	def check_infile(self, file, _format):
		_format = get_bioformat(file, _format)
		# try every available format
		fs = [_format] + self.formats if _format else self.formats
		content_file = None
		for _f in fs:
			try:
				content_file = self.read_infile(file, _f)
			except Exception as e:
				continue
			f = _f
			break
		if not content_file or not f:
			raise Exception()
		return content_file, f

	def prepare_export(self, outfile=None, _format=None) -> (str, str):
		outfile, f = super(BioService, self).prepare_export(outfile, _format)
		return outfile, get_bioformat(outfile, f)
