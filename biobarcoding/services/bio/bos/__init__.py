from .. import BioService
from ....rest import auth_filter
from ....db_models import DBSession


##
# CHADO BOS SERVICE
##

class BosService(BioService):

    def prepare_values(self, **values):

        values['sourceuri'] = values.get('sourceuri', values.get('filesAPI'))

        return super(BosService, self).prepare_values(**values)

    def pre_query(self, purpose='read'):
        from ....db_models.sysadmin import PermissionType
        from ....db_models.core import data_object_type_id
        from ....db_models.bioinformatics import FunctionalObject
        purpose_id = DBSession.query(PermissionType).filter(PermissionType.name == purpose).one().id
        bos_clause = DBSession.query(FunctionalObject.native_id) \
            .filter(auth_filter(FunctionalObject, purpose_id, [data_object_type_id[self.bos]]))
        bos_clause = [i for i, in bos_clause.all()]  # Cannot use .subquery(), because we have two -separate- databases

        from sqlalchemy import inspect
        bos_clause = inspect(self.orm).primary_key[0].in_(bos_clause)
        query = self.db.query(self.orm).filter(bos_clause)
        return query
