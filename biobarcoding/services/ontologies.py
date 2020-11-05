
def create_ontologies(name, definition = None, remote_url = None):
    return {'status':'success','message':'CREATE: ontology dummy completed.'}, 200

def read_ontologies(ontology_id = None, name = None):
    from biobarcoding.db_models import DBSessionChado
    from biobarcoding.db_models.chado import Cv
    result = DBSessionChado().query(Cv)
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
    onto_name = pronto.Ontology(input_file).metadata.default_namespace
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
            -d Pg -c {onto_name}'''])
    if err:
        import os
        return {'status':'failure','message':f'Ontology in {os.path.basename(input_file)} could not be imported.\n{err}'}, 500
    return {'status':'success','message':f'Ontology in {os.path.basename(input_file)} imported properly.\n{out}'}, 200

def export_ontologies(id):
    return {'status':'success','message':'EXPORT: ontology dummy completed'}, 200
