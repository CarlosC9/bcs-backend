def create_analysis(program, programversion, name = None, description = None, algorithm = None, sourcename = None, sourceversion = None, sourceuri = None, timeexecuted = None):
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
            timeexecuted = timeexecuted)
    return res


def read_analysis(analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if analysis_id:
        res = conn.analysis.get_analyses(analysis_id = analysis_id)
    else:
        res = conn.analysis.get_analyses()
    return res


def update_analysis(analysis_id, program, programversion, name = None, description = None, algorithm = None, sourcename = None, sourceversion = None, sourceuri = None, timeexecuted = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    # res = conn.analysis.update_analysis(
    #         program = program,
    #         programversion = programversion,
    #         name = name,
    #         description = description,
    #         algorithm = algorithm,
    #         sourcename = sourcename,
    #         sourceversion = sourceversion,
    #         sourceuri = sourceuri,
    #         timeexecuted = timeexecuted)
    return res

def delete_analysis(analysis_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if analysis_id:
        res = conn.analysis.delete_analyses(analysis_id = analysis_id)
    else:
        res = conn.analysis.delete_analyses()
    return res
