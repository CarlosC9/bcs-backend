from biobarcoding.db_models import DBSession as db_session
from biobarcoding.db_models import DBSessionChado as chado_session

from biobarcoding.rest import Issue, IType


def create_sequences(organism_id=None, analysis_id=None, residues=None):
    issues = [Issue(IType.WARNING, 'CREATE sequences: dummy completed')]
    content = { organism_id : organism_id, analysis_id : analysis_id, residues : residues }
    return issues, {k:v for k,v in content.items() if v is not None}, 200


def read_sequences(sequence_id=None, ids=None, organism_id=None, analysis_id=None, phylotree_id=None):
    content = { 'sequence_id':sequence_id, 'ids':' '.join(ids) if ids else None, 'organism_id':organism_id,
                'analysis_id':analysis_id, 'phylotree_id':phylotree_id }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        content = __get_query(sequence_id, ids, organism_id, analysis_id, phylotree_id)
        if sequence_id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ sequences: The sequences were successfully read')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ sequences: The sequences could not be read.')], 500
    return issues, content, status


def update_sequences(sequence_id, organism_id=None, analysis_id=None, residues=None):
    issues = [Issue(IType.WARNING, 'UPDATE sequences: dummy completed')]
    content = { sequence_id : sequence_id, organism_id : organism_id, analysis_id : analysis_id, residues : residues }
    return issues, {k:v for k,v in content.items() if v is not None}, 200


def __delete_from_bcs(feature_id):
    from biobarcoding.db_models.bioinformatics import Sequence
    db_session.query(Sequence).filter(Sequence.chado_feature_id == feature_id) \
        .delete(synchronize_session='fetch')


def delete_sequences(sequence_id=None, ids=None, organism_id=None, analysis_id=None):
    content = { sequence_id:sequence_id, 'ids':' '.join(ids) if ids else '', organism_id:organism_id, analysis_id:analysis_id }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        query = __get_query(sequence_id, ids, organism_id, analysis_id)
        for seq in query.all():
            __delete_from_bcs(seq.feature_id)
        resp = query.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE sequences: {resp} sequences were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'DELETE sequences: The sequences could not be removed.')], 500
    return issues, content, status


def import_sequences(input_file, organism_id=None, analysis_id=None, format='fasta'):
    content = { input_file : input_file, format : format, organism_id : organism_id, analysis_id : analysis_id }
    content = {k:v for k,v in content.items() if v is not None}
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if not organism_id:
        try:
            organism_id = conn.organism.add_organism(genus='organism',
                                                     species='undefined', common='', abbr='')['organism_id']
        except Exception as e:
            organism_id = conn.organism.get_organisms(species='undefined')[0]['organism_id']
    try:
        resp = conn.feature.load_fasta(input_file, organism_id, analysis_id=analysis_id, update=True)
        from Bio import SeqIO
        __seqs2bcs([seq.id for seq in SeqIO.parse(input_file, format)])
        issues, status = [Issue(IType.INFO, f'IMPORT sequences: {resp} sequences were successfully imported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'IMPORT sequences: file {input_file} could not be imported.')], 500
    return issues, content, status


def __seqs2bcs(names):
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


def export_sequences(sequence_id=None, ids=None, organism_id=None, analysis_id=None, output_file=None):
    if not output_file:
        output_file = "/tmp/output_seqs.fas"
    try:
        query = __get_query(sequence_id, ids, organism_id, analysis_id)
        with open(output_file, "w") as file:
            for seq in query.all():
                file.write(f'>{seq.uniquename}\n{seq.residues}\n')
        issues, status = [Issue(IType.INFO, f'EXPORT sequences: {query.count()} sequences were successfully exported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'EXPORT sequences: The sequences could not be exported.')], 500
    return issues, output_file, status


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
