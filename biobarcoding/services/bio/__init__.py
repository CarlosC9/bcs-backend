from ..main import BasicService


##
# CHADO TOOLS
##

class BioService(BasicService):

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

	def prepare_export(self, **kwargs):
		outfile, format = super(BioService, self).prepare_export(**kwargs)

		from .. import get_bioformat
		format = get_bioformat(outfile, format)

		return outfile, format
