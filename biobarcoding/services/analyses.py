def create_analyses(program, programversion, name = None, sourcename = None, description = None, algorithm = None, sourceversion = None, sourceuri = None, date_executed = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    try:
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
        return {'status':'success','message':f'The analysis "{program} {programversion}" created successfully.'}, 201
    except Exception as e:
        return {'status':'failure','message':f'The analysis "{program} {programversion}" could not be created.'}, 500

def read_analyses(analysis_id=None, name=None, program=None, programversion=None, algorithm=None, sourcename=None, sourceversion=None, sourceuri=None, description=None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    res = conn.analysis.get_analyses(
        analysis_id=analysis_id,
        name=name,
        program=program,
        programversion=programversion,
        algorithm=algorithm,
        sourcename=sourcename,
        sourceversion=sourceversion, sourceuri=sourceuri,
        description=description)
    if analysis_id:
        return res[0], 200
    return res, 200

def update_analyses(analysis_id, program, programversion, name = None, description = None, algorithm = None, sourcename = None, sourceversion = None, sourceuri = None, date_executed = None):
    return {'status':'success','message':'UPDATE: analysis dummy completed'}, 200

def delete_analyses(analysis_id, name=None, program=None, programversion=None, algorithm=None, sourcename=None, sourceversion=None, sourceuri=None, description=None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    try:
        res = conn.analysis.delete_analyses(
            analysis_id=analysis_id,
            name=name,
            program=program,
            programversion=programversion,
            algorithm=algorithm,
            sourcename=sourcename,
            sourceversion=sourceversion,
            sourceuri=sourceuri,
            description=description)
        return {'status':'success','message':f'{res} analyses were successfully removed.'}, 201
    except Exception as e:
        return {'status':'failure','message':f'The analyses could not be removed.'}, 500
