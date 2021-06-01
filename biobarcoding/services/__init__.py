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


def exec_cmds(*args):
    import subprocess
    out = err = []
    for cmd in args:
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


def get_or_create(session, model, **kwargs):
    instance = session.query(model).filter_by(**kwargs).first()
    if not instance:
        params = dict((k, v) for k, v in kwargs.items())
        instance = model(**params)
        session.add(instance)
        session.flush()
    return instance


# TODO in progress
def get_simple_query(session, model, **kwargs):
    kwargs = {k:v for k,v in kwargs.items() if v is not None}
    return session.query(model).filter_by(**kwargs)