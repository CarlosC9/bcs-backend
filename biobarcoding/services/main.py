from ..rest import Issue, IType
from . import log_exception, get_bioformat


def __getORM(entity):
    orm = None
    if entity == 'templates':
        from ..db_models.sysadmin import AnnotationFormTemplate as orm
    if entity == 'fields':
        from ..db_models.sysadmin import AnnotationFormField as orm
    if entity == 'template_fields':
        from ..db_models.sysadmin import AnnotationFormTemplateField as orm
    if entity == 'annotations':
        from ..db_models.sysadmin import AnnotationItem as orm
    return orm


def __getService(entity):
    service = None
    if entity == 'templates':
        from .annotation_forms.templates import AuxService as service
    if entity == 'fields':
        from .annotation_forms.fields import AuxService as service
    if entity == 'template_fields':
        from .annotation_forms.template_fields import AuxService as service
    if entity == 'annotations':
        from .annotation_forms.annotations import AuxService as service
    return service


def getCRUDIE(entity):

    ORM = __getORM(entity)
    AuxService = __getService(entity)

    class CRUDIE():

        issues = []
        content = None
        count = 0
        status = 200

        def create(self, **kwargs):
            try:
                values = AuxService.prepare_values(**kwargs)
                self.content = AuxService.create(**values)
                add = AuxService.add_create(self.content, **kwargs)
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


class AuxServiceSuper:

    def prepare_values(self):
        return None

    def create(self):
        return None

    def post_create(self):
        return None

    def read(self):
        return None

    def get_query(self):
        return None

    def aux_filter(self):
        return None

    def aux_order(self):
        return None

    def check_file(self):
        return None
