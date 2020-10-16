# phylotree.cvterm = 'taxonomy'
def create_taxonomies(name, comment = None):
    return {'status':'success','message':'dummy completed'}, 200

def read_taxonomies(id = None):
    from biobarcoding.db_models import DBSessionChado
    from biobarcoding.db_models.chado import Phylotree, Cvterm
    result = DBSessionChado().query(Phylotree)
    if id:
        result = result.filter(Phylotree.phylotree_id==id)
    response = []
    for value in result.all():
        tmp = value.__dict__
        tmp.pop('_sa_instance_state', None)
        response.append(tmp)
    return {'status':'success','message':response}, 200

def update_taxonomies(id, name = None, comment = None):
    return {'status':'success','message':'dummy completed'}, 200

def delete_taxonomies(id):
    return {'status':'success','message':'dummy completed'}, 200

def import_taxonomies(input_file, name = None, comment = None):
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

def export_taxonomies(id = None):
    return {'status':'success','message':'dummy completed'}, 200
