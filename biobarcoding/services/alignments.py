import os.path

from biobarcoding.db_models import DBSession as db_session
from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.rest import IType, Issue, filter_parse, paginator
from biobarcoding.services import get_or_create, log_exception
from biobarcoding.db_models.chado import Analysis, AnalysisCvterm, Cvterm, Organism, Feature, AnalysisFeature
from sqlalchemy import or_

from biobarcoding.services.analyses import __get_query_ordered, __aux_own_filter


def create(**kwargs):
    issues = [Issue(IType.WARNING, 'CREATE alignments: dummy completed')]
    return issues, None, 200


count = 0
def read(id=None, **kwargs):
    content = None
    try:
        content = __get_query(id, **kwargs)
        if id:
            # TODO: add taxa list and seqs len ?
            content = content.first()
        else:
            # TODO: add taxa list and seqs len ?
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ alignments: The alignments were successfully read.')], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, 'READ alignments: The alignments could not be read.')], 400
    return issues, content, count, status


def update(id, **kwargs):
    issues = [Issue(IType.WARNING, 'UPDATE alignments: dummy completed')]
    return issues, None, 200


def __delete_from_bcs(*args):
    from biobarcoding.db_models.bioinformatics import MultipleSequenceAlignment
    db_session.query(MultipleSequenceAlignment).filter(MultipleSequenceAlignment.chado_id.in_(*args)) \
        .delete(synchronize_session='fetch')


def delete(id=None, **kwargs):
    content = None
    try:
        # TODO: The BCS data are not being deleted yet.
        query = __get_query(id, **kwargs)
        from biobarcoding.services.sequences import delete as delete_sequences
        _ids = [msa.analysis_id for msa in query.all()]
        delete_sequences(filter=[{'analysis_id':{'op':'in','unary':_ids}}])
        __delete_from_bcs(_ids)
        resp = query.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE alignments: The {resp} alignments were successfully removed.')], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, 'DELETE alignments: The alignments could not be removed.')], 404
    return issues, content, status


def __seq_org(name):
    org = chado_session.query(Organism).join(Feature).filter(Feature.uniquename == name).first()
    if not org:
        org = get_or_create(chado_session, Organism, genus='organism', species='undefined')
        chado_session.add(org)
        chado_session.flush()
    return org


def __seq_cvterm(name):
    cvterm = chado_session.query(Cvterm).join(Feature).filter(Feature.uniquename == name).first()
    if not cvterm:
        cvterm = get_or_create(chado_session, Cvterm, name='contig')
    return cvterm


def __bind2src(feature, srcname):
    from biobarcoding.services.sequences import __get_query as get_sequences
    src = get_sequences(uniquename=srcname).first()
    if src:
        from biobarcoding.db_models.chado import Featureloc
        relationship = Featureloc(feature_id=feature.feature_id, srcfeature_id=src.feature_id)
        chado_session.add(relationship)


def __msafile2chado(input_file, msa, format):
    from Bio import AlignIO
    seqs = AlignIO.read(input_file, format or 'fasta')
    for seq in seqs:
        feature = Feature(uniquename=f'{seq.name}_msa{msa.analysis_id}', residues=f'{seq.seq}',
                          organism_id=__seq_org(seq.name).organism_id,
                          type_id=__seq_cvterm(seq.name).cvterm_id)
        chado_session.add(feature)
        chado_session.flush()
        __bind2src(feature, seq.name)
        chado_session.add(AnalysisFeature(analysis_id=msa.analysis_id, feature_id=feature.feature_id))
    return msa


def __msa2bcs(msa):
    from biobarcoding.db_models.bioinformatics import MultipleSequenceAlignment
    bcs_msa = get_or_create(db_session, MultipleSequenceAlignment, chado_id=msa.analysis_id)
    db_session.add(bcs_msa)
    db_session.flush()
    return bcs_msa


def import_file(input_file, format='fasta', **kwargs):
    content = None
    try:
        if not kwargs.get('program'):
            kwargs['program']='Multiple Sequence Alignment'
        if not kwargs.get('programversion'):
            kwargs['programversion']='(Imported file)'
        if not kwargs.get('name'):
            kwargs['name']=os.path.basename(input_file)
        msa = Analysis(**kwargs)
        chado_session.add(msa)
        chado_session.flush()
        __msa2bcs(msa)
        cvterm_id = chado_session.query(Cvterm.cvterm_id) \
            .filter(or_(Cvterm.name == 'Alignment', Cvterm.name == 'Sequence alignment')).first()
        msa_cvterm = get_or_create(chado_session, AnalysisCvterm, cvterm_id=cvterm_id, analysis_id=msa.analysis_id)
        chado_session.add(msa_cvterm)
        chado_session.flush()
        __msafile2chado(input_file, msa, format)
        issues, status = [Issue(IType.INFO,
                                f'IMPORT alignments: The {format} alignment were successfully imported.',
                                os.path.basename(input_file))], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR,
                                f'IMPORT alignments: The file {input_file} could not be imported.',
                                os.path.basename(input_file))], 409
    return issues, content, status


def export(id, format='fasta', value={}, **kwargs):
    content = None
    try:
        if format in ('fasta', 'nexus'):
            from biobarcoding.services.sequences import __get_query as get_seqs
            seqs = get_seqs(filter={'analysis_id': {'op': 'eq', 'unary': id}})
            from biobarcoding.services.sequences import __seqs2file as export_sequences
            content = export_sequences(seqs.all(), format=format, header_format=value.get('header'))
            issues, status = [Issue(IType.INFO, f'EXPORT alignments: The alignment was successfully imported.')], 200
        else:
            issues, status = [Issue(IType.ERROR, f'EXPORT alignments: The format {format} could not be exported.')], 404
        return issues, content, status
    except Exception as e:
        log_exception(e)
        return [Issue(IType.ERROR, f'EXPORT alignments: The alignment could not be exported.')], None, 404


def __get_query(id=None, **kwargs):
    msa_ids = chado_session.query(AnalysisCvterm.analysis_id) \
        .join(Cvterm, AnalysisCvterm.cvterm_id == Cvterm.cvterm_id) \
        .filter(or_(Cvterm.name == 'Alignment', Cvterm.name == 'Sequence alignment')).all()
    query = chado_session.query(Analysis).filter(Analysis.analysis_id.in_(msa_ids))
    global count
    count = 0
    if id:
        query = query.filter(Analysis.analysis_id == id)
    else:
        if 'filter' in kwargs:
            query = query.filter(filter_parse(Analysis, kwargs.get('filter'), __aux_own_filter))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            count = query.count()
            query = paginator(query, kwargs.get('pagination'))
    return query

# Alignment notation
def read_alignmentComments(self, id=None):
    return [Issue(IType.WARNING, 'READ sequences comment: dummy completed')], {}, 200


def create_alignmentComments(self):
    return [Issue(IType.WARNING, 'CREATE sequences comment: dummy completed')], {}, 200


def update_alignmentComments(self, id):
    return [Issue(IType.WARNING, 'UPDATE sequences comment: dummy completed')], {}, 200


def delete_alignmentComments(self, id):
    return [Issue(IType.WARNING, 'DELETE sequences comment: dummy completed')], {}, 200
