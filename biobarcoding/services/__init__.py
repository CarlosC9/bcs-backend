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


def orm2json(row):
    # TODO: missing inherited fields
    d = {}
    try:
        for column in row.__table__.columns:
            d[column.name] = str(getattr(row, column.name))
    except:
        pass
    return d


def get_or_create(session, model, **kwargs):
    instance = session.query(model).filter_by(**kwargs).first()
    if not instance:
        params = dict((k, v) for k, v in kwargs.items())
        instance = model(**params)
        session.add(instance)
        session.flush()
    return instance


def get_simple_query(session, model, **kwargs):
    kwargs = {k:v for k,v in kwargs.items() if v is not None}
    return session.query(model).filter_by(**kwargs)


def log_exception(e):
    import os, sys
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    print(e)