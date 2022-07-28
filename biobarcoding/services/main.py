import json
import os.path
import traceback

from flask import g
from sqlalchemy_continuum import transaction_class

from . import log_exception, get_orm_params, get_query, get_filtering
from ..rest import Issue, IType, filter_parse
from ..db_models import ORMBase, DBSession


def get_orm(entity):
    orm = None
    # SYS ORMS
    if entity == 'status_checkers':
        from ..db_models.jobs import StatusChecker as orm
    if entity == 'transactions':
        from ..db_models.chado import Transaction as orm
    # BOS ORMS
    elif entity == 'sequences':
        from ..db_models.chado import Feature as orm
    elif entity == 'alignments':
        from ..db_models.chado import Analysis as orm
    elif entity == 'phylotrees':
        from ..db_models.chado import Phylotree as orm
    elif entity == 'blasts':
        from ..db_models.bioinformatics import SequenceSimilarity as orm
    elif entity == 'discriminant_matrices':
        from ..db_models.bioinformatics import DiscriminantMatrix as orm
    elif entity == 'supermatrices':
        from ..db_models.bioinformatics import Supermatrix as orm
    elif entity == 'collections':
        from ..db_models.sysadmin import Collection as orm
    # BIO META ORMS
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
    # SYS SERVICES
    if entity == 'status_checkers':
        from .sys.status_checkers import Service
    if entity == 'transactions':
        from .sys.transactions import Service
    # BOS SERVICES
    elif entity == 'sequences':
        from .bio.bos.sequences import Service
    elif entity == 'alignments':
        from .bio.bos.alignments import Service
    elif entity == 'phylotrees':
        from .bio.bos.phylotrees import Service
    elif entity == 'blasts':
        from .bio.bos.blast import Service
    elif entity == 'discriminant_matrices':
        from .bio.bos.discriminant_matrices import Service
    elif entity == 'supermatrices':
        from .bio.bos.supermatrices import Service
    elif entity == 'analyses':
        from .bio.bos.analyses import Service
    elif entity == 'collections':
        from .bio.meta.collections import Service
    # BIO META SERVICES
    elif entity == 'individuals':
        from .bio.meta.individuals import Service
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
                traceback.print_exc()
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
        self.formats = []

    ##
    # CREATE
    ##

    # filter and deduce values (mostly ids) for creation or update
    def prepare_values(self, cv=None, cvterm=None, db=None, dbxref=None, **values):

        from ..db_models import DBSessionChado
        from ..db_models.chado import Cv, Cvterm, Db, Dbxref

        if cv and cvterm and not values.get('cvterm_id'):
            values['cvterm_id'] = DBSessionChado.query(Cvterm).join(Cv) \
                .filter(Cv.name == cv, Cvterm.name == cvterm).one().cvterm_id

        if db and dbxref and not values.get('dbxref_id'):
            values['dbxref_id'] = DBSessionChado.query(Dbxref).join(Db) \
                .filter(Db.name == db, Dbxref.accession == dbxref).one().dbxref_id

        if values.get('cvterm_id') and not values.get('dbxref_id'):
            values['dbxref_id'] = DBSessionChado.query(Cvterm) \
                .filter(Cvterm.cvterm_id == values.get('cvterm_id')) \
                .first().dbxref_id

        if values.get('dbxref_id') and not values.get('cvterm_id'):
            values['cvterm_id'] = DBSessionChado.query(Cvterm) \
                .filter(Cvterm.dbxref_id == values.get('dbxref_id')) \
                .first().cvterm_id

        return values

    # filter and deduce values (mostly ids) for creation or update
    def prepare_external_values(self, **values) -> dict:
        return values

    # check the validity of the values against the constraints to create or update
    def check_values(self, **values) -> dict:
        return get_orm_params(self.orm, **values)

    # deal with the creation
    def create(self, **kwargs) -> (object, int):
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
    def after_create(self, new_object, **values) -> dict:
        return values

    ##
    # READ
    ##

    # read like get_query, but telling if it asks only for one or more. It can use get_query
    # @return: result, total_count
    def read(self, **kwargs) -> (any, int):
        content, count = self.get_query(purpose='read', **kwargs)
        if kwargs.get('id'):
            try:
                content = self.attach_data(content.one())[0]
                return content, count
            except Exception as e:
                pass
        content = self.attach_data(*content.all())
        if kwargs.get('only_ids') or kwargs.get('values', {}).get('only_ids'):
            from sqlalchemy import inspect
            content = [c[inspect(self.orm).primary_key[0].name] for c in content]
        return content, count

    # any additional read if any
    def attach_data(self, *content) -> list:
        try:
            import json
            from ..common import generate_json
            return [json.loads(generate_json(c)) for c in content]      # TODO move to the place attachment happens
        except Exception as e:
            print('Error: The row could not be jsonify for attachment.')
            log_exception(e)
            return []

    ##
    # GET SQLALCHEMY QUERY
    ##

    # provide a sqlalchemy query (might be paged) and the total_count
    # @return: Query, total_count
    def get_query(self, query=None, id=None, purpose='delete', **kwargs) -> (object, int):
        _orm = self.orm
        if get_filtering('version_history', kwargs) or get_filtering('transaction_id', kwargs) \
                or get_filtering('version', kwargs) or get_filtering('version_id', kwargs) \
                or get_filtering('issued_at', kwargs):
            from sqlalchemy_continuum import version_class
            self.orm = version_class(self.orm)
        _resp = get_query(self.db, self.orm, query or self.pre_query(purpose), id,
                          aux_filter=self.aux_filter, aux_order=self.aux_order, **kwargs)
        self.orm = _orm
        return _resp

    # method to filter by acl and more particular issues when querying
    def pre_query(self, purpose) -> object:
        return None

    # method to filter by external values when querying
    def aux_filter(self, _filter: dict) -> list:
        clauses = []

        _v = _filter.get('version') or _filter.get('version_id')
        if _v:
            clauses.append(filter_parse(self.orm, {'transaction_id': _v}))

        _t = transaction_class(self.orm)
        if _filter.get("issued_at"):
            _ids = self.db.query(_t.id) \
                .filter(filter_parse(_t, {'issued_at': _filter.get("issued_at")}))
            clauses.append(self.orm.transaction_id.in_(_ids))

        if _filter.get("issued_at-from"):
            _ids = self.db.query(_t.id) \
                .filter(filter_parse(_t, {'issued_at': _filter.get("issued_at-from")}))
            clauses.append(self.orm.transaction_id.in_(_ids))

        if _filter.get("issued_at-to"):
            _ids = self.db.query(_t.id) \
                .filter(filter_parse(_t, {'issued_at': _filter.get("issued_at-to")}))
            clauses.append(self.orm.transaction_id.in_(_ids))

        return clauses

    # method to order by external values when querying
    def aux_order(self, order) -> list:
        return []

    ##
    # UPDATE
    ##

    # deal with the update
    def update(self, values={}, **kwargs) -> (any, int):
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
    def after_update(self, new_object, **values) -> dict:
        return values

    ##
    # DELETE
    ##

    # deal with the delete
    def delete(self, **kwargs) -> (any, int):
        content, count = self.get_query(purpose='delete', **kwargs)
        content = content.all()
        self.delete_related(*content, **kwargs)
        for row in content:
            self.db.delete(row)
        return content, count

    # any additional delete if any
    def delete_related(self, *content, **kwargs):
        return 0

    # TODO: generic import/export in progress
    ##
    # IMPORT
    ##

    # verify that a given file is the data it pretend to be
    def read_infile(self, file, _format) -> any:
        return file

    # figure out the format for a given file
    def check_infile(self, file, _format) -> (any, str):
        return self.read_infile(file, _format), _format

    def import_file(self, infile, format=None, **kwargs) -> (any, int):
        """
        content_file, _format = self.check_infile(infile, format)
        # values = self.prepare_values(**kwargs)
        content = self.create(**kwargs)
        c, count = self.file2data(content_file, **kwargs)
        return content, count
        """
        return None, 0

    ##
    # EXPORT
    ##

    def prepare_export(self, outfile=None, _format=None) -> (str, str):
        from flask import current_app
        from werkzeug.utils import secure_filename
        from biobarcoding import app_acronym
        outfile = os.path.join(current_app.config['UPLOAD_FOLDER'], secure_filename(outfile or f'output_{app_acronym}'))
        return outfile, _format

    def data2file(self, data: list, outfile, format: str, **kwargs) -> int:
        # TODO: check that the format is in self.formats ?
        with open(outfile, 'w') as wf:
            json.dump(data, wf)
        return len(data)

    def export_file(self, outfile=None, format=None, **kwargs) -> (str, int):
        outfile, _format = self.prepare_export(outfile, format)
        query, count = self.get_query(purpose='export', **kwargs)
        count = self.data2file(query.all(), outfile, _format, **kwargs)
        return outfile, count
