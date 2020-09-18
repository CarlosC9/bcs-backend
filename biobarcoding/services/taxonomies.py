# phylotree.cvterm = 'taxonomy'
def create_taxonomy(name, comment = None):
    return 'dummy completed'

def read_taxonomy(taxon_id = None):
    return 'dummy completed'

def update_taxonomy(taxon_id, name = None, comment = None):
    return 'dummy completed'

def delete_taxonomy(taxon_id = None):
    return 'dummy completed'

def import_taxonomy(input_file, organism_id = None):
    # yes '' | perl ./load_taxonomy_cvterms_edited.pl -H localhost -D postgres -u postgres -d Pg -p postgres;
    from flask import current_app
    with open(current_app.config["CHADO_CONF"], 'r') as chado_conf:
        import yaml
        # {cfg["host"]} {cfg["database"]} {cfg["user"]} {cfg["password"]} {cfg["schema"]} {cfg["port"]}
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
        return err
    return f'Taxa in {input_file} imported properly.'

def export_taxonomy(output_file = None, organism_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if not output_file:
        output_file = 'output_taxa.gbk'
    if not organism_id:
        organism_id = conn.organism.get_organisms(species='unknown')[0]['organism_id']
    import sys
    with open('/tmp/' + output_file, "w") as sys.stdout:
        resp = conn.export.export_gbk(organism_id)
        print(resp)
    return '/tmp/' + output_file
