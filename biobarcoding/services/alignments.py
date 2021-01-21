from biobarcoding.authentication import bcs_session
from biobarcoding.db_models import DBSession as db_session
from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.services import get_or_create
from biobarcoding.db_models.chado import Analysis, AnalysisCvterm, Cvterm, Organism, Feature, AnalysisFeature
from sqlalchemy import or_


@bcs_session(read_only=False)
def create_alignments(program=None, programversion=None, name=None, sourcename=None, description=None, algorithm=None,
                      sourceversion=None, sourceuri=None, timeexecuted=None):
    return {'status': 'success', 'message': f'CREATE: alignments dummy complete.'}, 201


@bcs_session(read_only=True)
def read_alignments(alignment_id=None, ids=None, name=None, program=None, programversion=None, algorithm=None,
                    sourcename=None, sourceversion=None, sourceuri=None, description=None, feature_id=None):
    result = __get_query(alignment_id=alignment_id, ids=ids, name=name, program=program, programversion=programversion,
                         algorithm=algorithm, sourcename=sourcename, sourceversion=sourceversion, sourceuri=sourceuri,
                         description=description, feature_id=feature_id)
    from biobarcoding.services import chado2json
    if alignment_id:
        return chado2json(result)[0], 200
    return chado2json(result), 200


@bcs_session(read_only=False)
def update_alignments(alignment_id, program, programversion, name=None, description=None, algorithm=None,
                      sourcename=None, sourceversion=None, sourceuri=None, timeexecuted=None):
    return {'status': 'success', 'message': 'UPDATE: alignments dummy completed'}, 200


def __delete_from_bcs(analysis_id):
    from biobarcoding.db_models.bioinformatics import MultipleSequenceAlignment
    db_session.query(MultipleSequenceAlignment).filter(MultipleSequenceAlignment.chado_analysis_id == analysis_id) \
        .delete(synchronize_session='fetch')


@bcs_session(read_only=False)
def delete_alignments(alignment_id=None, ids=None, name=None, program=None, programversion=None, algorithm=None,
                      sourcename=None, sourceversion=None, sourceuri=None, description=None):
    try:
        # TODO: The BCS data are not been deleted yet.
        res = __get_query(alignment_id=alignment_id, ids=ids, name=name, program=program, programversion=programversion,
                          algorithm=algorithm, sourcename=sourcename, sourceversion=sourceversion, sourceuri=sourceuri,
                          description=description)
        from biobarcoding.services.sequences import delete_sequences
        for msa in res.all():
            delete_sequences(analysis_id=msa.analysis_id)
            __delete_from_bcs(msa.analysis_id)
        res.delete(synchronize_session='fetch')
        return {'status': 'success', 'message': f'{res} alignments were successfully removed.'}, 201
    except Exception as e:
        print(e)
        return {'status': 'failure', 'message': f'The alignments could not be removed.'}, 500


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


def __msa2chado(input_file, msa, format):
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
    bcs_msa = get_or_create(db_session, MultipleSequenceAlignment, chado_analysis_id=msa.analysis_id)
    db_session.add(bcs_msa)
    db_session.flush()
    return bcs_msa


@bcs_session(read_only=False)
def import_alignments(input_file, format='fasta', **kwargs):
    msa = get_or_create(chado_session, Analysis, **kwargs)
    chado_session.add(msa)
    chado_session.flush()
    __msa2bcs(msa)
    cvterm_id = chado_session.query(Cvterm.cvterm_id) \
        .filter(or_(Cvterm.name == 'Alignment', Cvterm.name == 'Sequence alignment')).first()
    msa_cvterm = get_or_create(chado_session, AnalysisCvterm, cvterm_id=cvterm_id, analysis_id=msa.analysis_id)
    chado_session.add(msa_cvterm)
    chado_session.flush()
    try:
        __msa2chado(input_file, msa, format)
        return {'status': 'success', 'message': f'Imported MSA from {format} file.'}, 200
    except Exception as e:
        print(e)
        return {'status': 'failure', 'message': e}, 500


@bcs_session(read_only=True)
def export_alignments(analysis_id):
    from biobarcoding.services.sequences import export_sequences
    return export_sequences(analysis_id=analysis_id)


def __get_query(alignment_id=None, ids=None, name=None, program=None, programversion=None, algorithm=None,
                sourcename=None, sourceversion=None, sourceuri=None, description=None, feature_id=None):
    from biobarcoding.services.analyses import __get_query as get_analyses
    query = get_analyses(alignment_id, ids, name, program, programversion, algorithm,
                         sourcename, sourceversion, sourceuri, description, feature_id)
    msa_ids = chado_session.query(AnalysisCvterm.analysis_id) \
        .join(Cvterm, AnalysisCvterm.cvterm_id == Cvterm.cvterm_id) \
        .filter(or_(Cvterm.name == 'Alignment', Cvterm.name == 'Sequence alignment'))
    query = query.filter(Analysis.analysis_id.in_(msa_ids))
    return query


# Alignment notation
def read_alignmentComments(self, id=None):
    return {'status': 'success', 'message': f'Dummy completed'}, 200


def create_alignmentComments(self):
    return {'status': 'success', 'message': f'Dummy completed'}, 200


def update_alignmentComments(self, id):
    return {'status': 'success', 'message': f'Dummy completed'}, 200


def delete_alignmentComments(self, id):
    return {'status': 'success', 'message': f'Dummy completed'}, 200
