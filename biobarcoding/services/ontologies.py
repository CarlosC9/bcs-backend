from biobarcoding.db_models import DBSessionChado as chado_session


def create_ontologies(name, definition = None, remote_url = None):
    return {'status':'success','message':'CREATE: ontology dummy completed.'}, 200


def read_ontologies(ontology_id = None, name = None):
    from biobarcoding.db_models.chado import Cv
    result = chado_session.query(Cv)
    if ontology_id:
        result = result.filter(Cv.cv_id==ontology_id)
    if name:
        result = result.filter(Cv.name==name)
    response = []
    for value in result.all():
        tmp = value.__dict__
        tmp.pop('_sa_instance_state', None)
        response.append(tmp)
    if ontology_id:
        return response[0], 200
    return response, 200


def update_ontologies(ontology_id, name = None, definition = None, remote_url = None, input_file = None):
    return {'status':'success','message':'UPDATE: ontology dummy completed'}, 200


def delete_ontologies(ontology_id = None):
    return {'status':'success','message':'DELETE: ontology dummy completed'}, 200


def import_ontologies(input_file):
    from flask import current_app
    cfg = current_app.config
    # f"""go2fmt.pl -p obo_text -w xml {input_file} | \
    #     go-apply-xslt oboxml_to_chadoxml - > {input_file}.xml"""
    # cmd = f"""go2chadoxml {input_file} > /tmp/{input_file}.chado.xml;
    #     stag-storenode.pl -d 'dbi:Pg:dbname={cfg['database']};host={cfg['host']};port=cfg['port']'
    #     --user {cfg['user']} --password {cfg['password']} /tmp/{input_file}.chado.xml"""
    import pronto
    import os
    namespace = pronto.Ontology(input_file).metadata.default_namespace
    if namespace:
        onto_name = f'-c {namespace}'
    else:
        onto_name = f'-c {os.path.basename(input_file)}'
    from biobarcoding.services import exec_cmds
    out, err = exec_cmds([
        f'''perl ./biobarcoding/services/perl_scripts/gmod_load_cvterms.pl\
            -H {cfg["CHADO_HOST"]}\
            -D {cfg["CHADO_DATABASE"]}\
            -r {cfg["CHADO_USER"]}\
            -p {cfg["CHADO_PASSWORD"]}\
            -d Pg -s null -u\
            {input_file}''',
        f'''perl ./biobarcoding/services/perl_scripts/gmod_make_cvtermpath.pl\
            -H {cfg["CHADO_HOST"]}\
            -D {cfg["CHADO_DATABASE"]}\
            -u {cfg["CHADO_USER"]}\
            -p {cfg["CHADO_PASSWORD"]}\
            -d Pg {onto_name}'''])
    if err:
        return {'status':'failure','message':f'Ontology in {os.path.basename(input_file)} could not be imported.\n{err}'}, 500
    return {'status':'success','message':f'Ontology in {os.path.basename(input_file)} imported properly.\n{out}'}, 200


def export_ontologies(id):
    return {'status':'success','message':'EXPORT: ontology dummy completed'}, 200


# TODO:
#  featureprop, analysisprop, phylotreeprop ?
#  phylotree ?
def read_cvterms(cv_id = None, cvterm_id = None, feature_id=None, analysis_id=None, phylotree_id=None):
    from biobarcoding.db_models.chado import Cv, Cvterm
    result = chado_session.query(Cvterm)
    if cvterm_id:
        result = result.filter(Cvterm.cvterm_id==cvterm_id)
    if cv_id:
        result = result.filter(Cvterm.cv_id==cv_id)
    if feature_id:
        from biobarcoding.db_models.chado import FeatureCvterm
        cv_ids = chado_session.query(FeatureCvterm.cvterm_id)\
            .filter(FeatureCvterm.feature_id==feature_id)
        result = result.filter(Cvterm.cv_id.in_(cv_ids))
    if analysis_id:
        from biobarcoding.db_models.chado import AnalysisCvterm
        cv_ids = chado_session.query(AnalysisCvterm.cvterm_id)\
            .filter(AnalysisCvterm.analysis_id==analysis_id)
        result = result.filter(Cvterm.cv_id.in_(cv_ids))
    if phylotree_id:
        # from biobarcoding.db_models.chado import Phylotree
        # cv_id = chado_session.query(Phylotree.type_id)\
        #     .filter(Phylotree.phylotree_id==phylotree_id)
        # result = result.filter(Cvterm.cv_id==cv_id)
        pass
    response = []
    for value in result.all():
        tmp = value.__dict__
        tmp.pop('_sa_instance_state', None)
        response.append(tmp)
    if cvterm_id:
        return response[0], 200
    return response, 200
