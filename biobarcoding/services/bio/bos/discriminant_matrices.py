from . import BosService
from ...main import get_orm
from ....db_models import DBSession


##
# BLAST SERVICE
##
class Service(BosService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSession
        self.orm = get_orm('discriminant-matrices')
        self.bos = 'discriminant-matrix'