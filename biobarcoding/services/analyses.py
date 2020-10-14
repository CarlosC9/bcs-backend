def create_analysis(program, programversion, name = '', sourcename = '', description = None, algorithm = None, sourceversion = None, sourceuri = None, date_executed = None):
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
    return res, 0


def read_analysis(analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if analysis_id:
        res = conn.analysis.get_analyses(analysis_id = analysis_id)
    else:
        res = conn.analysis.get_analyses()
    return res, 0


def update_analysis(analysis_id, program, programversion, name = None, description = None, algorithm = None, sourcename = None, sourceversion = None, sourceuri = None, date_executed = None):
    return 'dummy completed', 0

def delete_analysis(analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if analysis_id:
        res = conn.analysis.delete_analyses(analysis_id = analysis_id)
    else:
        res = conn.analysis.delete_analyses()
    return res, 0
