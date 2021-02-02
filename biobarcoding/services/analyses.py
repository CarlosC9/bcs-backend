from biobarcoding.db_models import DBSessionChado as chado_session


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


def read_analyses(analysis_id=None, ids=None, name=None, program=None, programversion=None, algorithm=None, sourcename=None, sourceversion=None, sourceuri=None, description=None, feature_id=None):
    result = __get_query(
        analysis_id=analysis_id,
        ids=ids,
        name=name,
        program=program,
        programversion=programversion,
        algorithm=algorithm,
        sourcename=sourcename,
        sourceversion=sourceversion,
        sourceuri=sourceuri,
        description=description,
        feature_id=feature_id)
    from biobarcoding.services import chado2json
    if analysis_id:
        return chado2json(result)[0], 200
    return chado2json(result), 200


def update_analyses(analysis_id, program, programversion, name = None, description = None, algorithm = None, sourcename = None, sourceversion = None, sourceuri = None, date_executed = None):
    return {'status':'success','message':'UPDATE: analysis dummy completed'}, 200


def delete_analyses(analysis_id=None, ids=None, name=None, program=None, programversion=None, algorithm=None, sourcename=None, sourceversion=None, sourceuri=None, description=None):
    try:
        res = __get_query(
            analysis_id=analysis_id,
            ids=ids,
            name=name,
            program=program,
            programversion=programversion,
            algorithm=algorithm,
            sourcename=sourcename,
            sourceversion=sourceversion,
            sourceuri=sourceuri,
            description=description).delete(synchronize_session='fetch')
        return {'status':'success','message':f'{res} analyses were successfully removed.'}, 201
    except Exception as e:
        print(e)
        return {'status':'failure','message':f'The analyses could not be removed.'}, 500


def __get_query(analysis_id=None, ids=None, name=None, program=None, programversion=None, algorithm=None, sourcename=None, sourceversion=None, sourceuri=None, description=None, feature_id=None):
    from biobarcoding.db_models.chado import Analysis
    query = chado_session.query(Analysis)
    if analysis_id:
        query = query.filter(Analysis.analysis_id==analysis_id)
    if ids:
        query = query.filter(Analysis.analysis_id.in_(ids))
    if name:
        query = query.filter(Analysis.name==name)
    if program:
        query = query.filter(Analysis.program==program)
    if programversion:
        query = query.filter(Analysis.programversion==programversion)
    if algorithm:
        query = query.filter(Analysis.algorithm==algorithm)
    if sourcename:
        query = query.filter(Analysis.sourcename==sourcename)
    if sourceversion:
        query = query.filter(Analysis.sourceversion==sourceversion)
    if sourceuri:
        query = query.filter(Analysis.sourceuri==sourceuri)
    if description:
        query = query.filter(Analysis.description==description)
    if feature_id:
        from biobarcoding.db_models.chado import AnalysisFeature
        query = query.join(AnalysisFeature)\
            .filter(AnalysisFeature.feature_id==feature_id)
    return query
