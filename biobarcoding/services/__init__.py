def conn_chado():
    from flask import current_app
    from chado import ChadoInstance
    with open(current_app.config["CHADO_CONF"], 'r') as chado_conf:
        import yaml
        cfg = yaml.load(chado_conf, Loader=yaml.FullLoader)
        conn = ChadoInstance(dbhost=cfg["host"],
            dbname=cfg["database"],
            dbuser=cfg["user"],
            dbpass=cfg["password"],
            dbschema=cfg["schema"],
            dbport=cfg["port"])
    return conn


def exec_cmds(cmds):
    import subprocess
    out = err = []
    for cmd in cmds:
        print(cmd)
        process = subprocess.Popen(cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True)
        o, e = process.communicate()
        out.append(o)
        print(f'OUT: \n{o.decode("utf-8")}\n')
        if e:
            print(f'ERROR: \n{e.decode("utf-8")}')
            err.append(e)
    return out, err

# proc = subprocess.Popen('perl ./biobarcoding/services/perl_scripts/gmod_load_cvterms.pl -H localhost -D postgres -r postgres -p postgres -d Pg -s null -u /tmp/taxrank.obo',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
