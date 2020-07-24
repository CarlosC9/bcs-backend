import logging
import os
import sys

import redis
import sqlalchemy
from flask import Response
from sqlalchemy.pool import StaticPool

import biobarcoding
from biobarcoding.common import generate_json
from biobarcoding.common.gp_helpers import create_pg_database_engine, load_table
from biobarcoding.models import DBSession, ORMBase

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


def initialize_database_data():
    # Load base tables
    # load_table(DBSession, User, tm_default_users)
    # load_table(DBSession, Authenticator, tm_authenticators)
    # load_table(DBSession, CaseStudyStatus, tm_case_study_version_statuses)
    # load_table(DBSession, ObjectType, tm_object_types)
    # load_table(DBSession, PermissionType, tm_permissions)
    # # Create and insert a user
    # session = DBSession()
    # # Create test User, if it does not exist
    # u = session.query(User).filter(User.name == 'test_user').first()
    # if not u:
    #     u = User()
    #     u.name = "test_user"
    #     u.uuid = "27c6a285-dd80-44d3-9493-3e390092d301"
    #     session.add(u)
    #     session.commit()
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

