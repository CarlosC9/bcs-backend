# phylotree.cvterm = 'taxonomy'
def create_taxonomy(name, comment = None):
    return {'status':'success','message':'dummy completed'}, 200

def read_taxonomy(id = None):
    return {'status':'success','message':'dummy completed'}, 200

def update_taxonomy(id, name = None, comment = None):
    return {'status':'success','message':'dummy completed'}, 200

def delete_taxonomy(id):
    return {'status':'success','message':'dummy completed'}, 200

def import_taxonomy(input_file, name = None, comment = None):
    # yes '' | perl ./load_taxonomy_cvterms_edited.pl -H localhost -D postgres -u postgres -d Pg -p postgres;
    from flask import current_app
    with open(current_app.config["CHADO_CONF"], 'r') as chado_conf:
        import yaml
        cfg = yaml.load(chado_conf, Loader=yaml.FullLoader)
        cmd = f'(cd ./biobarcoding/services/taxonomy/ && perl ./load_ncbi_taxonomy.pl -H {cfg["host"]} -D {cfg["database"]} -u {cfg["user"]} -p {cfg["password"]} -d Pg -i {input_file})'
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
        return {'status':'failure','message':f'Taxonomy in {os.path.basename(input_file)} could not be imported.\n{err}'}, 500
    return {'status':'success','message':f'Taxonomy in {input_file} imported properly.'}, 200

def export_taxonomy(id = None):
    return {'status':'success','message':'dummy completed'}, 200
