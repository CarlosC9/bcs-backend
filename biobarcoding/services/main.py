from flask import g

from . import log_exception, get_bioformat, get_orm_params
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
    if entity == 'annotation_template':
        from ..db_models.sysadmin import AnnotationTemplate as orm
    if entity == 'annotation_field':
        from ..db_models.sysadmin import AnnotationField as orm
    if entity == 'annotation_text':
        from ..db_models.sysadmin import AnnotationText as orm
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
                g.commit_after = False
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
                g.commit_after = False
                self.issues += [Issue(IType.ERROR,
                                      f'READ {entity}: The {entity} could not be read.')]
                self.status = 400
            return self.issues, self.content, self.count, self.status

        def update(self, values={}, **kwargs):
            try:
                self.content = AuxService.update(values=values, **kwargs)
                extra = AuxService.after_update(self.content, values=values, **kwargs)
                self.issues += [Issue(IType.INFO,
                                      f'UPDATE {entity}: The {entity} was/were updated successfully.')]
                self.status = 200
            except Exception as e:
                log_exception(e)
                g.commit_after = False
                self.issues += [Issue(IType.ERROR,
                                      f'UPDATE {entity}: The {entity} could not be updated.')]
                self.status = 409
            return self.issues, self.content, self.count, self.status

        def delete(self, **kwargs):
            try:
                self.content, self.count = AuxService.delete(**kwargs)
                extra = AuxService.after_delete(self.content, **kwargs)
                self.issues += [Issue(IType.INFO,
                                      f'DELETE {entity}: The {entity} was/were removed successfully.')]
                self.status = 200
            except Exception as e:
                log_exception(e)
                g.commit_after = False
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
                g.commit_after = False
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
                g.commit_after = False
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
        from sqlalchemy.exc import SQLAlchemyError
        try:
            content = self.orm(**values)
            self.db.add(content)
            self.db.flush()
        except SQLAlchemyError as e:    # IntegrityError: already exists
            print(type(e))
            print(str(e.__dict__['orig']))
            raise e
        return content

    # any additional creation if any
    def after_create(self, new_object, **kwargs):
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
    def get_query(self, query=None, id=None, **kwargs):
        """
         reserved keywords in kwargs:
           'values': specific values of orm fields to filter
           'filter': advanced filtering clauses (see also filter_parse)
           'order': advanced ordering clauses (see also order_parse)
           'pagination': pageIndex and pageSize to paginate
           'searchValue': full-text search value (hopefully)
         otherwise it will be treated as 'value'
        """
        from ..rest import filter_parse, order_parse
        from ..services import paginator
        query = query or self.db.query(self.orm)
        count = 0
        if id:
            query = query.filter(self.orm.id == id)
        else:
            if not kwargs.get('values'):
                kwargs['values'] = {}
            for k, v in kwargs.items():
                if not k in ['values', 'filter', 'order', 'pagination', 'searchValue'] and v:
                    kwargs['values'][k] = v
            if kwargs.get('values'):
                query = query.filter_by(**self.prepare_values(**kwargs.get('values')))
            if kwargs.get('filter'):
                query = query.filter(filter_parse(self.orm, kwargs.get('filter'), self.aux_filter))
            if kwargs.get('order'):
                query = query.order_by(order_parse(self.orm, kwargs.get('order'), self.aux_order))
            count = query.count()
            if kwargs.get('pagination'):
                query = paginator(query, kwargs.get('pagination'))
        return query, count

    # method to filter by external values when querying
    def aux_filter(self, filter) -> list:
        return []

    # method to order by external values when querying
    def aux_order(self, order) -> list:
        return []

    # deal with the update
    def update(self, values={}, **kwargs):
        content, count = self.get_query(**kwargs)
        values = self.prepare_values(**values)
        # self.content = self.content.update(AuxService.prepare_values(**values))
        content = content.all()
        for row in content:
            for k, v in values.items():
                setattr(row, k, v)
        return content

    # any additional update if any
    def after_update(self, new_object, values={}, **kwargs):
        return None

    # deal with the delete
    def delete(self, **kwargs):
        content, count = self.get_query(**kwargs)
        # content = content.delete(synchronize_session='fetch')
        content = content.all()
        for row in content:
            DBSession.delete(row)
        return content, count

    # any additional delete if any
    def after_delete(self, new_object, **kwargs):
        return None

    # verify that a given file is the data it pretend to be
    def check_file(self, file) -> bool:
        return True
