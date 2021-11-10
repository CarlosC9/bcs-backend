from . import log_exception, get_bioformat, get_orm_params, get_query
from ..rest import Issue, IType
from ..db_models import ORMBase, DBSession


def get_orm(entity):
    orm = None
    if entity == 'templates':
        from ..db_models.sysadmin import AnnotationFormTemplate as orm
    if entity == 'fields':
        from ..db_models.sysadmin import AnnotationFormField as orm
    if entity == 'annotations':
        from ..db_models.sysadmin import AnnotationItem as orm
    return orm


def get_service(entity):
    AuxService = None
    if entity == 'templates':
        from .annotation_forms.templates import AuxService
    if entity == 'fields':
        from .annotation_forms.fields import AuxService
    if entity == 'annotations':
        from .annotation_forms.annotations import AuxService
    return AuxService()


def getCRUDIE(entity):

    ORM = get_orm(entity)
    AuxService = get_service(entity)

    class CRUDIE:

        issues = []
        content = None
        count = 0
        status = 200

        def __init__(self):
            pass

        def create(self, **kwargs):
            try:
                self.content = AuxService.create(**kwargs)
                extra = AuxService.after_create(self.content, **kwargs)
                self.issues += [Issue(IType.INFO,
                                      f'CREATE {entity}: The {entity} was created successfully.')]
                self.status = 201
            except Exception as e:
                log_exception(e)
                self.issues += [Issue(IType.ERROR,
                                      f'CREATE {entity}: The {entity} could not be created.')]
                self.status = 409
            return self.issues, self.content, self.count, self.status

        def read(self, id=None, **kwargs):
            try:
                self.content, self.count = AuxService.read(**kwargs)
                self.issues += [Issue(IType.INFO,
                                      f'READ {entity}: The {entity} were read successfully.')]
                self.status = 200
            except Exception as e:
                log_exception(e)
                self.issues += [Issue(IType.ERROR,
                                      f'READ {entity}: The {entity} could not be read.')]
                self.status = 400
            return self.issues, self.content, self.count, self.status

        def update(self, value={}, **kwargs):
            try:
                self.content, self.count = AuxService.get_query(**kwargs)
                self.content = self.content.update(AuxService.prepare_values(**value))
                self.issues += [Issue(IType.INFO,
                                      f'UPDATE {entity}: The {entity} was/were updated successfully.')]
                self.status = 200
            except Exception as e:
                log_exception(e)
                self.issues += [Issue(IType.ERROR,
                                      f'UPDATE {entity}: The {entity} could not be updated.')]
                self.status = 409
            return self.issues, self.content, self.count, self.status

        def delete(self, **kwargs):
            try:
                self.content, self.count = AuxService.get_query(**kwargs)
                self.content = self.content.delete(synchronize_session='fetch')
                self.issues += [Issue(IType.INFO,
                                      f'DELETE {entity}: The {entity} was/were removed successfully.')]
                self.status = 200
            except Exception as e:
                log_exception(e)
                self.issues += [Issue(IType.ERROR,
                                      f'DELETE {entity}: The {entity} could not be removed.')]
                self.status = 404
            return self.issues, self.content, self.count, self.status

        # TODO: generic import in progress
        def import_file(self, input_file, format=None, **kwargs):
            format = get_bioformat(input_file, format)
            try:
                # AuxService.check_file(input_file, format)
                values = AuxService.prepare_values(**kwargs)
                # create required rows ?
                # file2db ?
                self.issues += [Issue(IType.INFO,
                                      f'IMPORT {entity}: The file {input_file} was imported successfully.')]
                self.status = 200
            except Exception as e:
                log_exception(e)
                self.issues += [Issue(IType.ERROR,
                                      f'IMPORT {entity}: The file {input_file} could not be imported.')]
                self.status = 409
            return self.issues, self.content, self.count, self.status

        # TODO: generic export in progress
        def export_file(self, format='fasta', **kwargs):
            try:
                # get_query ?
                # export ?
                self.issues = [Issue(IType.INFO,
                                        f'EXPORT {entity}: The {entity} were exported successfully.')]
                self.status = 200
            except Exception as e:
                log_exception(e)
                self.issues += [Issue(IType.ERROR,
                                      f'EXPORT {entity}: The {entity} could not be exported.')]
                self.status = 409
            return self.issues, self.content, self.count, self.status

    return CRUDIE()


class SimpleAuxService:

    def __init__(self):
        self.orm = ORMBase
        self.db = DBSession

    # filter and deduce values (mostly ids) for creation or update
    def prepare_values(self, **kwargs) -> dict:
        return get_orm_params(self.orm, **kwargs)

    # check the validity of the values against the constraints to create or update
    def check_values(self, **kwargs) -> dict:
        return kwargs

    # deal with the creation
    def create(self, **kwargs):
        values = self.prepare_values(**kwargs)
        content = self.orm(**values)
        self.db.add(content)
        self.db.flush()
        return content

    # any additional creation if any
    def after_create(self, new_obj, **kwargs):
        return None

    # read like get_query, but telling if it asks only for one or more. It can use get_query
    # @return: result, total_count
    def read(self, **kwargs):
        content, count = self.get_query(**kwargs)
        if kwargs.get('id'):
            content = content.first()
        else:
            content = content.all()
        return content, count

    # provide a sqlalchemy query (might be paged) and the total_count
    # @return: Query, total_count
    def get_query(self, **kwargs):
        return get_query(self.db, self.orm,
                         aux_filter=self.aux_filter,
                         aux_order=self.aux_order, **kwargs)

    # method to filter by external values when querying
    def aux_filter(self, filter) -> list:
        return []

    # method to order by external values when querying
    def aux_order(self, order) -> list:
        return []

    # verify that a given file is the data it pretend to be
    def check_file(self, file) -> bool:
        return True
