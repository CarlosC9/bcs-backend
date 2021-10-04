import os.path

from Bio import AlignIO

from ..db_models import DBSession as db_session
from ..db_models import DBSessionChado as chado_session
from ..db_models.chado import Organism, Feature, AnalysisFeature, Analysis
from ..db_models.bioinformatics import MultipleSequenceAlignment

from ..rest import IType, Issue
from ..services import get_or_create, log_exception, orm2json, get_bioformat
from ..services.analyses import __get_query as get_ansis_query
from ..services.sequences import __get_query as get_seqs_query, \
    create as create_seq, \
    delete as delete_sequences, \
    __seqs2file as export_sequences


##
# CREATE
##

def __msa2bcs(msa):
    return get_or_create(db_session, MultipleSequenceAlignment,
                         chado_id=msa.analysis_id,
                         chado_table='analysis',
                         name=msa.name)


def create(**kwargs):
    content = None
    try:
        if not kwargs.get('program'):
            kwargs['program'] = 'Multiple sequence alignment'
        # create as analysis
        from .analyses import create as create_ansis
        issues, content, status = create_ansis(**kwargs)
        chado_session.flush()   # analysis_id required
        msa_bcs = __msa2bcs(content)
        if status < 300:
            issues += [Issue(IType.INFO, f'CREATE alignments: The alignment "{kwargs.get("name")}" was created successfully.')]
        else:
            raise Exception(*issues)
    except Exception as e:
        log_exception(e)
        issues += [Issue(IType.ERROR, f'CREATE alignments: The alignment "{kwargs.get("name")}" could not be created.')]
        status = 409
    return issues, content, status


##
# READ
##

def read(id=None, **kwargs):
    content, count = None, 0
    try:
        content, count = __get_query(id, **kwargs)
        if id:
            content = orm2json(content.one())
            from sqlalchemy.sql.expression import func, distinct
            info = chado_session.query(func.max(func.length(Feature.residues)),
                                       func.count(AnalysisFeature.analysis_id),
                                       func.array_agg(distinct(Organism.genus + ' ' + Organism.species))) \
                .select_from(AnalysisFeature).join(Feature).join(Organism) \
                .filter(AnalysisFeature.analysis_id == content['analysis_id'])\
                .group_by(AnalysisFeature.analysis_id)
            content['seqlen'], content['seqnum'], content['taxa'] = info.one()
        else:
            # TODO: add taxa list, seqs len and seqs num ?
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ alignments: The alignments were successfully read.')], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, 'READ alignments: The alignments could not be read.')], 400
    return issues, content, count, status


##
# UPDATE
##

def update(id, **kwargs):
    from .analyses import update as update_ansis
    return update_ansis(id, **kwargs)


##
# DELETE
##

def __delete_from_bcs(*ids):
    query = db_session.query(MultipleSequenceAlignment).filter(MultipleSequenceAlignment.chado_id.in_(ids))
    return query.count()
    # TODO: check why all rows are deleted in bcs without filtering
    # return query.delete(synchronize_session='fetch')


def delete(id=None, **kwargs):
    content = None
    try:
        # TODO: The BCS data are not being deleted yet.
        content, count = __get_query(id, **kwargs)
        ids = [msa.analysis_id for msa in content.all()]
        delete_sequences(filter=[{'analysis_id':{'op':'in','unary':ids}}])
        bcs_delete = __delete_from_bcs(*ids)
        content = content.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE alignments: The {content} alignments were successfully removed.')], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, 'DELETE alignments: The alignments could not be removed.')], 404
    return issues, content, status


##
# IMPORT
##

def __seq_org_id(name):
    try:
        return chado_session.query(Organism).join(Feature).filter(Feature.uniquename == name).one().organism_id
    except Exception as e:
        return None


def __bind2src(feature, srcname):
    try:
        src = get_seqs_query(uniquename=srcname)[0].one()
        from ..db_models.chado import Featureloc
        relationship = Featureloc(feature_id=feature.feature_id, srcfeature_id=src.feature_id)
        chado_session.add(relationship)
    except Exception as e:
        relationship = None
    return relationship


def __msafile2chado(msa, seqs):
    for seq in seqs:
        issues, feature, status = create_seq(
            uniquename=f'{seq.id}.msa{msa.analysis_id}',
            name=seq.id,
            stock=seq.id,
            residues=str(seq.seq),
            organism_id=__seq_org_id(seq.id),
            type='sequence', subtype='aligned')
        __bind2src(feature, seq.id)
        chado_session.add(AnalysisFeature(analysis_id=msa.analysis_id, feature_id=feature.feature_id))
    return msa


def import_file(input_file, format=None, **kwargs):
    content = None
    format = get_bioformat(input_file, format)
    try:
        # check aligned file
        content_file = AlignIO.read(input_file, format or 'fasta')
        # Set missing default values
        if not kwargs.get('programversion'):
            kwargs['programversion'] = '(Imported file)'
        if not kwargs.get('sourcename'):
            kwargs['sourcename'] = os.path.basename(input_file)
        # Analysis row could exist for jobs, so get or create
        try:
            unique_keys = ['job_id'] if kwargs.get('job_id') else ('program', 'programversion', 'sourcename')
            content, count = get_ansis_query(**{k:kwargs[k] for k in unique_keys if k in kwargs})
            content = content.one()
            # Analysis row could be created by importing the results of other jobs, and without register as msa in bcs
            msa_bcs = __msa2bcs(content)
        except Exception as e:
            issues, content, status = create(**kwargs)
        # Read and import file content
        __msafile2chado(content, content_file)
        issues, status = [Issue(IType.INFO,
                                f'IMPORT alignments: The {format} alignment were successfully imported.',
                                os.path.basename(input_file))], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR,
                                f'IMPORT alignments: The file {input_file} could not be imported.',
                                os.path.basename(input_file))], 409
    return issues, content, status


##
# EXPORT
##

def export(id, format='fasta', value={}, **kwargs):
    content = None
    try:
        if format in ('fasta', 'nexus'):
            seqs = get_seqs_query(filter={'analysis_id': {'op': 'eq', 'unary': id}})[0]
            content = export_sequences(seqs.all(), format=format, header_format=value.get('header'), only_headers=value.get('only_headers'))
            issues, status = [Issue(IType.INFO, f'EXPORT alignments: The alignment was successfully imported.')], 200
        else:
            issues, status = [Issue(IType.ERROR, f'EXPORT alignments: The format {format} could not be exported.')], 404
        return issues, content, status
    except Exception as e:
        log_exception(e)
        return [Issue(IType.ERROR, f'EXPORT alignments: The alignment could not be exported.')], None, 404


##
# GETTER AND OTHERS
##

def __get_query(id=None, **kwargs):
    aln_clause = db_session.query(MultipleSequenceAlignment.chado_id).all()
    aln_clause = [i for i, in aln_clause]
    aln_clause = Analysis.analysis_id.in_(aln_clause)
    query = chado_session.query(Analysis).filter(aln_clause)
    return get_ansis_query(id, **kwargs, query=query)


##
# Alignment notation
##

def read_alignmentComments(self, id=None):
    return [Issue(IType.WARNING, 'READ alignments comment: dummy completed')], {}, 200


def create_alignmentComments(self):
    return [Issue(IType.WARNING, 'CREATE alignments comment: dummy completed')], {}, 200


def update_alignmentComments(self, id):
    return [Issue(IType.WARNING, 'UPDATE alignments comment: dummy completed')], {}, 200


def delete_alignmentComments(self, id):
    return [Issue(IType.WARNING, 'DELETE alignments comment: dummy completed')], {}, 200
