from .. import BioService
from ....rest import auth_filter
from ....db_models import DBSession


##
# CHADO BOS SERVICE
##

class BosService(BioService):

    def pre_query(self, purpose='read'):
        from ....db_models.sysadmin import PermissionType
        from ....db_models.core import data_object_type_id
        from ....db_models.bioinformatics import FunctionalObject
        purpose_id = DBSession.query(PermissionType).filter(PermissionType.name == purpose).one().id
        seq_clause = DBSession.query(FunctionalObject.native_id) \
            .filter(auth_filter(FunctionalObject, purpose_id, [data_object_type_id[self.bos]]))
        seq_clause = [i for i, in seq_clause.all()]

        from sqlalchemy import inspect
        seq_clause = inspect(self.orm).primary_key[0].in_(seq_clause)
        query = self.db.query(self.orm).filter(seq_clause)
        return query
