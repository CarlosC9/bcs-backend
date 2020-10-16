def create_analyses(program, programversion, name = '', sourcename = '', description = None, algorithm = None, sourceversion = None, sourceuri = None, date_executed = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    res = conn.analysis.add_analysis(
            program = program,
            programversion = programversion,
            name = name,
            description = description,
            algorithm = algorithm,
            sourcename = sourcename,
            sourceversion = sourceversion,
            sourceuri = sourceuri,
            date_executed = date_executed)
    return {'status':'success','message':res}, 0

def read_analyses(analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if analysis_id:
        res = conn.analysis.get_analyses(analysis_id = analysis_id)
    else:
        res = conn.analysis.get_analyses()
    return {'status':'success','message':res}, 0


def update_analyses(analysis_id, program, programversion, name = None, description = None, algorithm = None, sourcename = None, sourceversion = None, sourceuri = None, date_executed = None):
    return {'status':'success','message':'dummy completed'}, 0

def delete_analyses(analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if analysis_id:
        res = conn.analysis.delete_analyses(analysis_id = analysis_id)
    else:
        res = conn.analysis.delete_analyses()
    return {'status':'success','message':res}, 0
