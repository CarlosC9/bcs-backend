from .. import get_bioformat
from ..main import BasicService
from ...db_models import DBSessionChado


##
# CHADO TOOLS
##

class BioService(BasicService):

	def __init__(self):
		super(BioService, self).__init__()
		self.db = DBSessionChado

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
