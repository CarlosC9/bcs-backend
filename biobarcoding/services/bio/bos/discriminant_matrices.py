from . import BosService
from ...main import get_orm
from ....db_models import DBSession
from ....db_models.bioinformatics import DiscriminantMatrix


##
# BLAST SERVICE
##
class Service(BosService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSession
        self.orm = get_orm('discriminant_matrices')
        self.obj_type = 'discriminant-matrix'
        self.fos = DiscriminantMatrix
