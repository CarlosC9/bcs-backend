def create_sequences(organism_id = None, analysis_id = None, residues = None):
    return {'status':'success','message':'CREATE: sequences dummy completed'}, 200


def read_sequences(sequence_id = None, organism_id = None, analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    resp = conn.feature.get_features(organism_id = organism_id,
        analysis_id = analysis_id, name = sequence_id)
    if sequence_id:
        return resp[0], 200
    return resp, 200


def update_sequences(sequence_id, organism_id = None, analysis_id = None, residues = None):
    return {'status':'success','message':'UPDATE: sequences dummy completed'}, 200


def delete_sequences(sequence_id = None, organism_id = None, analysis_id = None, ids = None):
    resp = __get_query(sequence_id, organism_id, analysis_id, ids).delete(synchronize_session='fetch')
    from biobarcoding.db_models import DBSessionChado
    DBSessionChado.commit()
    return {'status':'success','message':f'{resp} sequences were successfully removed.'}, 200


def import_sequences(input_file, organism_id = None, analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if not organism_id:
        organism_id = conn.organism.get_organisms(species='unknown')[0]['organism_id']
    try:
        # resp = conn.feature.load_fasta(input_file, organism_id, analysis_id=analysis_id, sequence_type='polypeptide', update=True)
        resp = conn.feature.load_fasta(input_file, organism_id, analysis_id=analysis_id, update=True)
        return {'status':'success','message':f'Sequences: {resp}'}, 200
    except Exception as e:
        return {'status':'failure','message':e}, 500


def export_sequences(sequence_id = None, organism_id = None, analysis_id = None, ids = None, output_file = None):
    if not output_file:
        output_file = "output_seqs.fas"
    query = __get_query(sequence_id, organism_id, analysis_id, ids)
    import sys
    with open('/tmp/' + output_file, "w") as file:
        for seq in query.all():
            file.write(f'>{seq.uniquename}\n{seq.residues}\n')
    return '/tmp/' + output_file, 200

def __get_query(sequence_id = None, organism_id = None, analysis_id = None, ids = None):
    from biobarcoding.db_models import DBSessionChado
    from biobarcoding.db_models.chado import Feature
    query = DBSessionChado.query(Feature)
    if sequence_id:
        query = query.filter(Feature.feature_id==sequence_id)
    if organism_id:
        query = query.filter(Feature.organism_id==organism_id)
    if analysis_id:
        query = query.filter(Feature.analysis_id==analysis_id)
    if ids:
        query = query.filter(Feature.feature_id.in_(ids))
    return query
