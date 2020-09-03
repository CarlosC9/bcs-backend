def create_sequence(organism_id = None, analysis_id = None, residues = None):
    return 'dummy completed'

def get_sequence(sequence_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if sequence_id:
        resp = conn.feature.get_features(sequence_id)
    else:
        resp = conn.feature.get_features()
    return resp

def update_sequence(sequence_id, organism_id = None, analysis_id = None, residues = None):
    return 'dummy completed'

def delete_sequence(sequence_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if sequence_id:
        resp = conn.feature.delete_features(sequence_id)
    else:
        resp = conn.feature.delete_features()
    return resp

def import_sequences(input_file, organism_id = None, analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if not organism_id:
        organism_id = conn.organism.get_organisms(species='unknown')[0]['organism_id']
    if not analysis_id:
        analysis_id = conn.analysis.get_analyses(name='Unknown analysis')[0]['analysis_id']
    try:
        resp = conn.feature.load_fasta(input_file, organism_id, update=True)
    except Exception as e:
        resp = e
    return resp

def export_sequences(output_file = None, organism_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if not output_file:
        output_file = "output_seqs.fas"
    if not organism_id:
        organism_id = conn.organism.get_organisms(species='unknown')[0]['organism_id']
    import sys
    with open('/tmp/' + output_file, "w") as sys.stdout:
        conn.export.export_fasta(organism_id)
    return '/tmp/' + output_file
