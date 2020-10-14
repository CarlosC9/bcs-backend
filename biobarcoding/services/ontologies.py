def create_ontologies(name, definition = None, remote_url = None):
    return {'status':'success','message':'dummy completed'}, 200

def read_ontologies(ontology_id = None, name = None):
    return {'status':'success','message':'dummy completed'}, 200

def update_ontologies(ontology_id, name = None, definition = None, remote_url = None, input_file = None):
    return {'status':'success','message':'dummy completed'}, 200

def delete_ontologies(ontology_id = None):
    return {'status':'success','message':'dummy completed'}, 200

def import_ontologies(input_file, name = None, definition = None):
    from flask import current_app
    with open(current_app.config["CHADO_CONF"], 'r') as chado_conf:
        import yaml
        cfg = yaml.load(chado_conf, Loader=yaml.FullLoader)
        # f"""go2fmt.pl -p obo_text -w xml {input_file} | \
        #     go-apply-xslt oboxml_to_chadoxml - > {input_file}.xml"""
        cmd = f"""go2chadoxml {input_file} > /tmp/{input_file}.chado.xml;
            stag-storenode.pl -d 'dbi:Pg:dbname={cfg['database']};host={cfg['host']};port=cfg['port']'
            --user {cfg['user']} --password {cfg['password']} /tmp/{input_file}.chado.xml"""
    import subprocess
    process = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         shell=True)
    out, err = process.communicate()
    print(f'OUT: {out}\n')
    if err:
        print(err)
        import os
        return {'status':'failure','message':f'Ontology in {os.path.basename(input_file)} could not be imported.\n{err}'}, 500
    return {'status':'success','message':f'Ontology in {input_file} imported properly.'}, 200

def export_ontologies(id = None, name = None):
    return {'status':'success','message':'dummy completed'}, 200
