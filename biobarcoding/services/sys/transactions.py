from . import SysService
from ..main import get_orm
from ...db_models import DBSessionChado


##
# TRANSACTIONS TOOLS
##

class Service(SysService):

	def __init__(self):
		super(Service, self).__init__()
		self.db = DBSessionChado
		self.orm = get_orm('transactions')
