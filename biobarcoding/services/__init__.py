def conn_chado():
    from flask import current_app
    from chado import ChadoInstance
    conn = ChadoInstance(dbhost=current_app.config["CHADO_HOST"],
        dbname=current_app.config["CHADO_DATABASE"],
        dbuser=current_app.config["CHADO_USER"],
        dbpass=current_app.config["CHADO_PASSWORD"],
        dbschema=current_app.config["CHADO_SCHEMA"],
        dbport=current_app.config["CHADO_PORT"])
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


def chado2json(query):
    response = []
    for value in query.all():
        tmp = value.__dict__
        tmp.pop('_sa_instance_state', None)
        response.append(tmp)
    return response
