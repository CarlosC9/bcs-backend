from .analyses import Service as AnsisService
from ...main import get_orm
from ....db_models import DBSession
from ....db_models.bioinformatics import SequenceSimilarity


##
# BLAST SERVICE
##
class Service(AnsisService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSession
        self.orm = get_orm('blasts')
        self.obj_type = 'sequence-similarity'
        self.fos = SequenceSimilarity
