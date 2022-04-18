import json
import os.path

from flask import g
from typing import Tuple

from . import log_exception, get_bioformat, get_orm_params
from ..rest import Issue, IType
from ..db_models import ORMBase, DBSession


def get_orm(entity):
    orm = None
    # BOS ORMS
    # TODO: add discriminant-matrices, blasts, supermatrices, collections
    if entity == 'sequences':
        from ..db_models.chado import Feature as orm
    elif entity == 'alignments':
        from ..db_models.chado import Analysis as orm
    elif entity == 'phylotrees':
        from ..db_models.chado import Phylotree as orm
    # BIO META ORMS
    # TODO: add publications, topics, sources, crs (reports?)
    elif entity == 'individuals':
        from ..db_models.chado import Stock as orm
    elif entity == 'analyses':
        from ..db_models.chado import Analysis as orm
    elif entity == 'taxonomies':
        from ..db_models.chado import Phylotree as orm
    elif entity == 'organisms':
        from ..db_models.chado import Organism as orm
    elif entity == 'ontologies':
        from ..db_models.chado import Cv as orm
    elif entity == 'cvterms':
        from ..db_models.chado import Cvterm as orm
    # ANNOTATION ORMS
    elif entity in ('template', 'templates'):
        from ..db_models.sa_annotations import AnnotationFormTemplate as orm
    elif entity in ('field', 'fields'):
        from ..db_models.sa_annotations import AnnotationFormField as orm
    elif entity in ('annotation', 'annotations'):
        from ..db_models.sa_annotations import AnnotationItem as orm
    elif entity == 'annotation_template':
        from ..db_models.sa_annotations import AnnotationTemplate as orm
    elif entity == 'annotation_field':
        from ..db_models.sa_annotations import AnnotationField as orm
    elif entity == 'annotation_text':
        from ..db_models.sa_annotations import AnnotationText as orm
    return orm


def get_service(entity):
    Service = None
    # BOS SERVICES
    # TODO: add discriminant-matrices, blasts, supermatrices, collections
    if entity == 'sequences':
        from .bio.bos.sequences import Service
    elif entity == 'alignments':
        from .bio.bos.alignments import Service
    elif entity == 'phylotrees':
        from .bio.bos.phylotrees import Service
    # BIO META SERVICES
    # TODO: add publications, topics, sources, crs (reports?)
    elif entity == 'individuals':
        from .bio.meta.individuals import Service
    elif entity == 'analyses':
        from .bio.meta.analyses import Service
    elif entity == 'taxonomies':
        from .bio.meta.taxonomies import Service
    elif entity == 'organisms':
        from .bio.meta.organisms import Service
    elif entity == 'ontologies':
        from .bio.meta.ontologies import Service
    elif entity == 'cvterms':
        from .bio.meta.ontologies import CvtermService as Service
    # ANNOTATION SERVICES
    elif entity == 'templates':
        from .annotation_forms.templates import Service
    elif entity == 'fields':
        from .annotation_forms.fields import Service
    elif entity == 'annotations':
        from .annotation_forms.annotations import Service
    return Service()


def getCRUDIE(entity):

    ORM = get_orm(entity)
    Service = get_service(entity)

    class CRUDIE:

        issues = []
        content = None
        count = 0
        status = 200

        def __init__(self):
            pass

        def create(self, **kwargs):
            try:
                self.content, self.count = Service.create(**kwargs)
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

        def read(self, **kwargs):
            try:
                self.content, self.count = Service.read(**kwargs)
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

        def update(self, **kwargs):
            try:
                self.content, self.count = Service.update(**kwargs)
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
                self.content, self.count = Service.delete(**kwargs)
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

        def import_data(self, input_file, **kwargs):
            try:
                self.content, self.count = Service.import_file(input_file, **kwargs)
                self.issues += [Issue(IType.INFO,
                                      f'IMPORT {entity}: The file {os.path.basename(input_file)} was imported successfully.')]
                self.status = 200
            except Exception as e:
                log_exception(e)
                g.commit_after = False
                self.issues += [Issue(IType.ERROR,
                                      f'IMPORT {entity}: The file {os.path.basename(input_file)} could not be imported.')]
                self.status = 409
            return self.issues, self.content, self.count, self.status

        def export_data(self, **kwargs):
            try:
                self.content, self.count = Service.export_file(**kwargs)
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


class BasicService:

    def __init__(self):
        self.orm = ORMBase
        self.db = DBSession

    ##
    # CREATE
    ##

    # filter and deduce values (mostly ids) for creation or update
    def prepare_values(self, **values) -> dict:
        return values

    # filter and deduce values (mostly ids) for creation or update
    def prepare_external_values(self, **values) -> dict:
        return values

    # check the validity of the values against the constraints to create or update
    def check_values(self, **values) -> dict:
        return get_orm_params(self.orm, **values)

    # deal with the creation
    def create(self, **kwargs) -> Tuple[any, int]:
        values = self.prepare_values(**kwargs)
        from sqlalchemy.exc import SQLAlchemyError
        try:
            content = self.orm(**self.check_values(**values))
            self.db.add(content)
            self.db.flush()
        except SQLAlchemyError as e:    # IntegrityError: already exists
            print(type(e))
            print(str(e.__dict__.get('orig', f'Exception "{e}" (unknown orig)')))
            raise e
        foreign_values = self.prepare_external_values(**kwargs)
        self.after_create(content, **foreign_values)
        return content, 1

    # any additional creation if any
    def after_create(self, new_object, **values):
        return new_object

    ##
    # READ
    ##

    # read like get_query, but telling if it asks only for one or more. It can use get_query
    # @return: result, total_count
    def read(self, **kwargs) -> Tuple[any, int]:
        content, count = self.get_query(purpose='read', **kwargs)
        if kwargs.get('id'):
            content = self.attach_data(content.first())
        else:
            content = [self.attach_data(c) for c in content.all()]
        return content, count

    # any additional read if any
    def attach_data(self, content):
        try:
            import json
            from ..common import generate_json
            return json.loads(generate_json(content))
        except:
            return None

    ##
    # GET SQLALCHEMY QUERY
    ##

    # provide a sqlalchemy query (might be paged) and the total_count
    # @return: Query, total_count
    def get_query(self, query=None, id=None, purpose='delete', **kwargs) -> Tuple[object, int]:
        """
         reserved keywords in kwargs:
           'values': specific values of orm fields to filter
           'filter': advanced filtering clauses (see also filter_parse)
           'order': advanced ordering clauses (see also order_parse)
           'pagination': pageIndex and pageSize to paginate
           'searchValue': full-text search value (hopefully)
         otherwise it will be treated as 'values'
        """
        from ..rest import filter_parse, order_parse
        from ..services import paginator
        query = query or self.pre_query(purpose) or self.db.query(self.orm)
        count = 0
        if id:
            if hasattr(self.orm, 'id'):
                query = query.filter(self.orm.id == id)
            else:
                from sqlalchemy import inspect
                query = query.filter(inspect(self.orm).primary_key[0] == id)
            count = query.count()
        else:
            if not kwargs.get('values'):
                kwargs['values'] = {}
            for k, v in kwargs.items():
                if not k in ['values', 'filter', 'order', 'pagination', 'searchValue'] and v:
                    kwargs['values'][k] = v
            if kwargs.get('values'):
                query = query.filter_by(**get_orm_params(self.orm, **kwargs.get('values')))
                # query = query.filter(filter_parse(self.orm, kwargs.get('values'), self.aux_filter))
            if kwargs.get('filter'):
                query = query.filter(filter_parse(self.orm, kwargs.get('filter'), self.aux_filter))
            if kwargs.get('searchValue', '') != '':
                if hasattr(self.orm, "ts_vector"):
                    query = query.filter(self.orm.ts_vector.match(kwargs.get('searchValue')))
            if kwargs.get('order'):
                query = query.order_by(*order_parse(self.orm, kwargs.get('order'), self.aux_order))
            count = query.count()
            if kwargs.get('pagination'):
                query = paginator(query, kwargs.get('pagination'))
        return query, count

    # method to filter by acl and more particular issues when querying
    def pre_query(self, purpose):
        return None

    # method to filter by external values when querying
    def aux_filter(self, filter) -> list:
        return []

    # method to order by external values when querying
    def aux_order(self, order) -> list:
        return []

    ##
    # UPDATE
    ##

    # deal with the update
    def update(self, values={}, **kwargs) -> Tuple[any, int]:
        changes = self.prepare_values(**values)
        content, count = self.get_query(purpose='contribute', **kwargs)
        foreign_changes = self.prepare_external_values(**values)
        content = content.all()
        for row in content:
            for k, v in changes.items():
                try:
                    setattr(row, k, v)
                except:
                    log_exception(f'"{k}" does not exist in "{row}"')
            self.after_update(row, **foreign_changes)
        return content, count

    # any additional update if any
    def after_update(self, new_object, **values):
        return new_object

    ##
    # DELETE
    ##

    # deal with the delete
    def delete(self, **kwargs) -> Tuple[any, int]:
        content, count = self.get_query(purpose='delete', **kwargs)
        # content = content.delete(synchronize_session='fetch')
        content = content.all()
        for row in content:
            self.db.delete(row)
        self.after_delete(*content, **kwargs)
        return content, count

    # any additional delete if any
    def after_delete(self, *content, **kwargs):
        # TODO: check why all rows are deleted in bcs without filtering
        return content

    # TODO: generic import/export in progress
    ##
    # IMPORT
    ##

    # verify that a given file is the data it pretend to be
    def check_file(self, file, format) -> bool:
        return True

    def import_file(self, infile, format=None, **kwargs) -> Tuple[any, int]:
        """
        format = get_bioformat(file, format)
        self.check_file(file, format)
        values = self.prepare_values(**kwargs)
        """
        # create required rows ?
        # file2db ?
        return None, 0

    ##
    # EXPORT
    ##

    def prepare_export(self, outfile=None, format=None, **kwargs) -> Tuple[str, str]:
        from flask import current_app
        from werkzeug.utils import secure_filename
        from biobarcoding import app_acronym
        outfile = os.path.join(current_app.config['UPLOAD_FOLDER'], secure_filename(outfile or f'output_{app_acronym}'))
        return outfile, format

    def data2file(self, data: list, outfile, format: str, **kwargs) -> int:
        # TODO: check that the format is in self.formats ?
        with open(outfile, 'w') as wf:
            json.dump(data, wf)
        return len(data)

    def export_file(self, outfile=None, format=None, **kwargs) -> Tuple[any, int]:
        outfile, format = self.prepare_export(outfile=outfile, format=format)
        query, count = self.get_query(purpose='export', **kwargs)
        count = self.data2file(query.all(), outfile, format, **kwargs)
        return outfile, count
