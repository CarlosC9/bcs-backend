def create_sequences(organism_id = None, analysis_id = None, residues = None):
    return {'status':'success','message':'dummy completed'}, 200

def read_sequences(sequence_id = None, organism_id = None, analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    resp = conn.feature.get_features(organism_id = organism_id,
        analysis_id = analysis_id, name = sequence_id)
    return {'status':'success','message':resp}, 200

def update_sequences(sequence_id, organism_id = None, analysis_id = None, residues = None):
    return {'status':'success','message':'dummy completed'}, 200

def delete_sequences(sequence_id = None, organism_id = None, analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    resp = conn.feature.delete_features(uniquename = sequence_id, organism_id = organism_id, analysis_id = analysis_id)
    return {'status':'success','message':resp}, 200

def import_sequences(input_file, organism_id = None, analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if not organism_id:
        organism_id = conn.organism.get_organisms(species='unknown')[0]['organism_id']
    if not analysis_id:
        analysis_id = conn.analysis.get_analyses(name='Unknown analysis')[0]['analysis_id']
    try:
        resp = conn.feature.load_fasta(input_file, organism_id, analysis_id=analysis_id, update=True)
        return {'status':'success','message':resp}, 200
    except Exception as e:
        return {'status':'failure','message':e}, 500

def export_sequences(output_file = None, sequence_id = None, organism_id = None, analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if not output_file:
        output_file = "output_seqs.fas"
    if not organism_id:
        organism_id = conn.organism.get_organisms(species='unknown')[0]['organism_id']
    import sys
    stdout = sys.stdout
    with open('/tmp/' + output_file, "w") as sys.stdout:
        conn.export.export_fasta(organism_id)
    sys.stdout = stdout
    return '/tmp/' + output_file, 200