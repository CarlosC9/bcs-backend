from biobarcoding.authentication import bcs_session
from biobarcoding.db_models import DBSession as db_session
from biobarcoding.db_models import DBSessionChado as chado_session


@bcs_session(read_only=False)
def create_sequences(organism_id=None, analysis_id=None, residues=None):
    return {'status': 'success', 'message': 'CREATE: sequences dummy completed'}, 200


@bcs_session(read_only=True)
def read_sequences(sequence_id=None, ids=None, organism_id=None, analysis_id=None, phylotree_id=None):
    # from biobarcoding.services import conn_chado
    # conn = conn_chado()
    # resp = conn.feature.get_features(organism_id = organism_id,
    #     analysis_id = analysis_id, name = sequence_id)
    from biobarcoding.services import chado2json
    if sequence_id:
        return chado2json(__get_query(sequence_id, ids, organism_id, analysis_id, phylotree_id))[0], 200
    return chado2json(__get_query(sequence_id, ids, organism_id, analysis_id, phylotree_id)), 200


@bcs_session(read_only=False)
def update_sequences(sequence_id, organism_id=None, analysis_id=None, residues=None):
    return {'status': 'success', 'message': 'UPDATE: sequences dummy completed'}, 200


def __delete_from_bcs(feature_id):
    from biobarcoding.db_models.bioinformatics import Sequence
    db_session.query(Sequence).filter(Sequence.chado_feature_id == feature_id) \
        .delete(synchronize_session='fetch')


@bcs_session(read_only=False)
def delete_sequences(sequence_id=None, ids=None, organism_id=None, analysis_id=None):
    query = __get_query(sequence_id, ids, organism_id, analysis_id)
    for seq in query.all():
        __delete_from_bcs(seq.feature_id)
    resp = query.delete(synchronize_session='fetch')
    return {'status': 'success', 'message': f'{resp} sequences were successfully removed.'}, 200


@bcs_session(read_only=False)
def import_sequences(input_file, organism_id=None, analysis_id=None, format='fasta'):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if not organism_id:
        try:
            organism_id = conn.organism.add_organism(genus='organism',
                                                     species='undefined', common='', abbr='')['organism_id']
        except Exception as e:
            organism_id = conn.organism.get_organisms(species='undefined')[0]['organism_id']
    try:
        # resp = conn.feature.load_fasta(input_file, organism_id, analysis_id=analysis_id, sequence_type='polypeptide', update=True)
        resp = conn.feature.load_fasta(input_file, organism_id, analysis_id=analysis_id, update=True)
        from Bio import SeqIO
        __ids2bcs([seq.id for seq in SeqIO.parse(input_file, format)])
        return {'status': 'success', 'message': f'Sequences: {resp}'}, 200
    except Exception as e:
        print(e)
        return {'status': 'failure', 'message': e}, 500


def __ids2bcs(names):
    from biobarcoding.db_models.chado import Feature
    from biobarcoding.db_models import DBSession as db_session
    for seq in chado_session.query(Feature).filter(Feature.uniquename.in_(names)).all():
        db_session.merge(__feature2bcs(seq))


def __feature2bcs(seq):
    from biobarcoding.services import get_or_create
    from biobarcoding.db_models.bioinformatics import Specimen, Sequence
    bcs_specimen = get_or_create(db_session, Specimen, name=seq.uniquename)
    db_session.merge(bcs_specimen)
    db_session.flush()
    bcs_sequence = get_or_create(db_session, Sequence,
                                 chado_feature_id=seq.feature_id,
                                 chado_table='feature',
                                 name=seq.uniquename,
                                 specimen_id=bcs_specimen.id)
    return bcs_sequence


@bcs_session(read_only=True)
def export_sequences(sequence_id=None, ids=None, organism_id=None, analysis_id=None, output_file=None):
    if not output_file:
        output_file = "output_seqs.fas"
    query = __get_query(sequence_id, ids, organism_id, analysis_id)
    import sys
    with open('/tmp/' + output_file, "w") as file:
        for seq in query.all():
            file.write(f'>{seq.uniquename}\n{seq.residues}\n')
    return '/tmp/' + output_file, 200


def __get_query(sequence_id=None, ids=None, organism_id=None, analysis_id=None, phylotree_id=None, uniquename=None):
    from biobarcoding.db_models.chado import Feature
    query = chado_session.query(Feature)
    if sequence_id:
        query = query.filter(Feature.feature_id == sequence_id)
    if ids:
        query = query.filter(Feature.feature_id.in_(ids))
    if organism_id:
        query = query.filter(Feature.organism_id == organism_id)
    if uniquename:
        query = query.filter(Feature.uniquename == uniquename)
    if analysis_id:
        from biobarcoding.db_models.chado import AnalysisFeature
        analysis_ids = chado_session.query(AnalysisFeature.feature_id).filter(AnalysisFeature.analysis_id==analysis_id).all()
        query = query.filter(Feature.feature_id.in_(analysis_ids))
    if phylotree_id:
        from biobarcoding.db_models.chado import Phylonode
        phylotree_ids = chado_session.query(Phylonode.feature_id).filter(Phylonode.phylotree_id==phylotree_id).all()
        query = query.filter(Feature.feature_id.in_(phylotree_ids))
    return query
