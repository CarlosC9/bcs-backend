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
