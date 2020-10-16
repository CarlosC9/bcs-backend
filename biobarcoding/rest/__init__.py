import logging
import os
import sys

import redis
import sqlalchemy
from flask import Response, Blueprint
from sqlalchemy import orm
from sqlalchemy.pool import StaticPool

import biobarcoding
from biobarcoding.common import generate_json
from biobarcoding.common.gp_helpers import create_pg_database_engine, load_table
from biobarcoding.db_models import DBSession, ORMBase, DBSessionChado, ORMBaseChado, ObjectType
from biobarcoding.db_models.bioinformatics import *
from biobarcoding.db_models.sysadmin import *
from biobarcoding.db_models.geographics import *
from biobarcoding.db_models.metadata import *
from biobarcoding.db_models.jobs import *

bcs_api_base = "/api"  # Base for all RESTful calls
bcs_gui_base = "/gui"  # Base for the Angular2 GUI
bcs_external_gui_base = "/gui_external"  # Base for the Angular2 GUI when loaded from URL

logger = logging.getLogger(__name__)
logging.getLogger('flask_cors').level = logging.DEBUG
log_level = logging.DEBUG


def build_json_response(obj, status=200):
    return Response(generate_json(obj),
                    mimetype="text/json",
                    status=status)


def prepare_default_configuration(create_directories):
    def default_directories(path, tmp_path):
        REDIS_HOST = "redis"
        REDIS_PORT = 6379
        BROKER_URL = f"redis://{REDIS_HOST}:{REDIS_PORT}/0"
        BACKEND_URL = BROKER_URL
        return {
            "DB_CONNECTION_STRING": f"sqlite:///{path}/bcs.db",
            "CACHE_FILE_LOCATION": f"{tmp_path}/cache",
            "CELERY_BROKER_URL": BROKER_URL,
            "CELERY_BACKEND_URL": BACKEND_URL,
            "REDIS_HOST_FILESYSTEM_DIR": f"{tmp_path}/sessions",
            "GOOGLE_APPLICATION_CREDENTIALS": f"{path}/firebase-key.json",
            "CHADO_CONF": f"{path}/chado_conf.yml"
        }

    from appdirs import AppDirs
    app_dirs = AppDirs("bcs-backend")

    # Default directories, multi-platform
    data_path = app_dirs.user_data_dir
    cache_path = app_dirs.user_cache_dir
    if create_directories:
        os.makedirs(data_path, exist_ok=True)
        os.makedirs(cache_path, exist_ok=True)

    # Obtain and create directories
    dirs = default_directories(data_path, cache_path)
    for v in dirs.values():
        if v:
            if create_directories:
                os.makedirs(v, exist_ok=True)

    # Default configuration
    return f"""{os.linesep.join([f'{k}="{v}"' for k, v in dirs.items()])}
# Flask Session (server side session)
REDIS_HOST="filesystem:local_session"
TESTING="True"
SELF_SCHEMA=""
""", data_path + os.sep + "bcs_local.conf"


def load_configuration_file(flask_app):
    # Initialize configuration
    config_file_var = "BCS_CONFIG_FILE"
    try:
        _, file_name = prepare_default_configuration(False)
        if not os.environ.get(config_file_var):
            found = False
            for f in [file_name]:
                if os.path.isfile(f):
                    print(f"Assuming {f} as configuration file")
                    found = True
                    os.environ[config_file_var] = f
                    break
            if not found:
                cfg, file_name = prepare_default_configuration(True)
                print(f"Generating {file_name} as configuration file:\n{cfg}")
                with open(file_name, "wt") as f:
                    f.write(cfg)
                os.environ[config_file_var] = file_name

        print("-----------------------------------------------")
        print(f'Configuration file at: {os.environ[config_file_var]}')
        print("-----------------------------------------------")
        flask_app.config.from_envvar(config_file_var)
    except Exception as e:
        print(f"{config_file_var} environment variable not defined!")
        print(e)


def construct_session_persistence_backend(flask_app):
    # A REDIS instance needs to be available. Check it
    # A local REDIS could be as simple as:
    #
    # docker run --rm -p 6379:6379 redis:alpine
    #
    d = {}
    if 'REDIS_HOST' in flask_app.config:
        r_host = flask_app.config['REDIS_HOST']
        d["SESSION_KEY_PREFIX"] = "nis:"
        d["SESSION_PERMANENT"] = False
        rs2 = None
        if r_host == "redis_lite":
            try:
                import redislite
                rs2 = redislite.Redis("tmp_bcs_backend_redislite.db")  # serverconfig={'port': '6379'}
                d["SESSION_TYPE"] = "redis"
                d["SESSION_REDIS"] = rs2
                # d["PERMANENT_SESSION_LIFETIME"] = 3600
            except ImportError as e:
                print("Package 'redislite' not found. Please, either change REDIS_HOST configuration variable to "
                      "'filesystem' or 'redis', or execute 'pip install redislite' and retry")
                sys.exit(1)
        elif r_host.startswith("filesystem:"):
            d["SESSION_TYPE"] = "filesystem"
            if flask_app.config.get("REDIS_HOST_FILESYSTEM_DIR"):
                d["SESSION_FILE_DIR"] = flask_app.config.get("REDIS_HOST_FILESYSTEM_DIR")
            d["SESSION_FILE_THRESHOLD"] = 100
            # d["SESSION_FILE_MODE"] = 666
        else:
            rs2 = redis.Redis(r_host)
            d["SESSION_TYPE"] = "redis"
            d["SESSION_REDIS"] = rs2
            # d["PERMANENT_SESSION_LIFETIME"] = 3600
        if rs2:
            try:
                print("Trying connection to REDIS '"+r_host+"'")
                rs2.ping()
                print("Connected to REDIS instance '"+r_host+"'")
            except:
                print("REDIS instance '"+r_host+"' not reachable, exiting now!")
                sys.exit(1)
        elif "SESSION_TYPE" not in d:
            print("No session persistence backend configured, exiting now!")
            sys.exit(1)
    return d


# "22bc2577-883a-4408-bba9-fcace20c0fc8":
# "e80a7d27-3ec8-4aa1-b49c-5498e0f85bee":
# "d30120f0-28df-4bca-90e4-4f0676d1c874":
# "83084df6-7ad0-45d7-b3f1-6de594c78611":
# "7e23991b-24a0-4da1-8251-c3c3434dfb87":
# "bfc0c9fe-631f-44d0-8e96-e22d01ffb1ed":
# "dec7e901-b3f4-4343-b3d1-4fa5fbf3804e":
# "013b2f3b-4b2f-4b6c-8f5f-425132aea74b":
# "3eef41be-fde1-4ad4-92d0-fe795158b41d":
# "0fba3591-4ffc-4a88-977a-6e1d922f0735":
# "a61fc587-1272-4d46-bdd0-027cde1b8a78":
# "600397ef-0102-486e-a6f7-d43b0f8ce4b9":
# "763e57b7-2636-4c04-9861-d865fe0bb5ab":
# "788065a7-d9f5-46fa-b8ba-8bc223d09331":
# "38fb34f7-a952-4036-9b0b-4d6c59e8f8d4":
# "0292821a-dd33-450a-bdd8-813b2b95c456":

tm_object_types = {  # ObjectType
    "22bc2577-883a-4408-bba9-fcace20c0fc8": "sequence",
    "e80a7d27-3ec8-4aa1-b49c-5498e0f85bee": "multiple-sequence-alignment",
    "d30120f0-28df-4bca-90e4-4f0676d1c874": "phylogenetic-tree",
    "83084df6-7ad0-45d7-b3f1-6de594c78611": "geographic-layer",
    "7e23991b-24a0-4da1-8251-c3c3434dfb87": "sequence-alignment"
}

tm_permissions = {  # PermissionType
    "f19ad19f-0a74-44e8-bd4e-4762404a35aa": "read",
    "04cac7ca-a90b-4d12-a966-d8b0c77fca70": "annotate",
    "d0924822-32fa-4456-8143-0fd48da33fd7": "contribute",
    "83d837ab-01b2-4260-821b-8c4a3c52e9ab": "share",
    "d3137471-84a0-4bcf-8dd8-16387ea46a30": "delete"
}

tm_default_users = {  # Identities
    "f3848599-4aa3-4964-b7e1-415d478560be": "admin",
    "2a512320-cef7-41c6-a141-8380d900761b": "_anonymous",
    "27c6a285-dd80-44d3-9493-3e390092d301": "test_user",
}

tm_authenticators = {  # Authenticator
    "b33193c3-63b9-49f7-b888-ceba619d2812": "firebase-google",
    "c09fa36b-62a3-4904-9600-e3bb5028d809": "firebase-facebook",
    "f510cb30-7a44-4cb1-86f5-1b112e43293a": "firebase-mail",
    "5f32a593-306f-4b69-983c-0a5680556fae": "local",
}


def initialize_database_data():
    # Load base tables
    load_table(DBSession, Identity, tm_default_users)
    load_table(DBSession, Authenticator, tm_authenticators)
    load_table(DBSession, ObjectType, tm_object_types)
    load_table(DBSession, PermissionType, tm_permissions)
    DBSession.remove()


def initialize_database(flask_app):
    recreate_db = False
    if 'DB_CONNECTION_STRING' in flask_app.config:
        db_connection_string = flask_app.config['DB_CONNECTION_STRING']
        print("Connecting to BCS database server")
        print(db_connection_string)
        print("-----------------------------")
        if db_connection_string.startswith("sqlite://"):
            biobarcoding.engine = sqlalchemy.create_engine(db_connection_string,
                                                         echo=True,
                                                         connect_args={'check_same_thread': False},
                                                         poolclass=StaticPool)
        else:
            biobarcoding.engine = create_pg_database_engine(db_connection_string, "bcs", recreate_db=recreate_db)

        # global DBSession # global DBSession registry to get the scoped_session
        DBSession.configure(bind=biobarcoding.engine)  # reconfigure the sessionmaker used by this scoped_session
        orm.configure_mappers()  # Important for SQLAlchemy-Continuum
        tables = ORMBase.metadata.tables
        connection = biobarcoding.engine.connect()
        table_existence = [biobarcoding.engine.dialect.has_table(connection, tables[t].name) for t in tables]
        connection.close()
        if False in table_existence:
            ORMBase.metadata.bind = biobarcoding.engine
            ORMBase.metadata.create_all()
        # connection = biobarcoding.engine.connect()
        # table_existence = [biobarcoding.engine.dialect.has_table(connection, tables[t].name) for t in tables]
        # connection.close()

        # Load base tables
        initialize_database_data()
    else:
        print("No database connection defined (DB_CONNECTION_STRING), exiting now!")
        sys.exit(1)


def initialize_database_chado(flask_app):
    if 'CHADO_CONF' in flask_app.config:
        with open(flask_app.config["CHADO_CONF"], 'r') as chado_conf:
            import yaml
            cfg = yaml.load(chado_conf, Loader=yaml.FullLoader)
        db_connection_string = f'postgres://{cfg["user"]}:{cfg["password"]}@{cfg["host"]}:{cfg["port"]}/{cfg["database"]}'
        print("Connecting to Chado database server")
        print(db_connection_string)
        print("-----------------------------")
        biobarcoding.chado_engine = sqlalchemy.create_engine(db_connection_string, echo=True)
        # global DBSessionChado # global DBSessionChado registry to get the scoped_session
        DBSessionChado.configure(bind=biobarcoding.chado_engine)  # reconfigure the sessionmaker used by this scoped_session
        orm.configure_mappers()  # Important for SQLAlchemy-Continuum
        ORMBaseChado.metadata.bind = biobarcoding.chado_engine
        ORMBaseChado.metadata.reflect()
    else:
        print("No database connection defined (DB_CONNECTION_STRING), exiting now!")
        sys.exit(1)


def register_api(bp: Blueprint, view, endpoint: str, url: str, pk='id', pk_type='int'):
    view_func = view.as_view(endpoint)
    bp.add_url_rule(url, defaults={pk: None}, view_func=view_func, methods=['GET'])
    bp.add_url_rule(url, view_func=view_func, methods=['POST'])
    bp.add_url_rule(f'{url}<{pk_type}:{pk}>', view_func=view_func, methods=['GET', 'PUT', 'DELETE'])
