from ...main import get_orm
from ....db_models import DBSession
from .alignments import Service as AlgnService


##
# BLAST SERVICE
##
class Service(AlgnService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSession
        self.orm = get_orm('supermatrices')
        self.bos = 'supermatrix'