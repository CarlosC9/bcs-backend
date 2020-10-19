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
from biobarcoding.db_models import DBSession, ORMBase, ObjectType
from biobarcoding.db_models.bioinformatics import *
from biobarcoding.db_models.sysadmin import *
from biobarcoding.db_models.geographics import *
from biobarcoding.db_models.metadata import *
from biobarcoding.db_models.jobs import *
from biobarcoding.db_models.hierarchies import *

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
            "GOOGLE_APPLICATION_CREDENTIALS": f"{path}/firebase-key.json"
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

tm_job_statuses = {
    "22bc2577-883a-4408-bba9-fcace20c0fc8": "created",
    "e80a7d27-3ec8-4aa1-b49c-5498e0f85bee": "preparing_workspace",
    "763e57b7-2636-4c04-9861-d865fe0bb5ab": "exporting_to_supported_file_formats",
    "d30120f0-28df-4bca-90e4-4f0676d1c874": "transferring_data_to_resource",
    "83084df6-7ad0-45d7-b3f1-6de594c78611": "submitting",
    "a61fc587-1272-4d46-bdd0-027cde1b8a78": "waiting_for_execution",
    "600397ef-0102-486e-a6f7-d43b0f8ce4b9": "executing",
    "bfc0c9fe-631f-44d0-8e96-e22d01ffb1ed": "transferring_data_from_resource",
    "788065a7-d9f5-46fa-b8ba-8bc223d09331": "importing_into_database",
    "7e23991b-24a0-4da1-8251-c3c3434dfb87": "completed_successfully",
    "0fba3591-4ffc-4a88-977a-6e1d922f0735": "completed_error",
    "dec7e901-b3f4-4343-b3d1-4fa5fbf3804e": "cleaning_up_workspace",
    "013b2f3b-4b2f-4b6c-8f5f-425132aea74b": "cancelled",
    "3eef41be-fde1-4ad4-92d0-fe795158b41d": "paused",
}

tm_job_mgmt_types = {
    "38fb34f7-a952-4036-9b0b-4d6c59e8f8d4": "galaxy",
    "0292821a-dd33-450a-bdd8-813b2b95c456": "ssh"
}

tm_processes = {
    "02f44e54-f139-4ea0-a1bf-fe27054c0d6c": "klustal-1",
    "903a73a9-5a4e-4cec-b8fa-4fc9bd5ffab5": "blast",
    "5c4ba6db-e7f2-4d5c-a89a-76059ac116b1": "mrbayes",
    "ea647c4e-2063-4246-bd9a-42f6a57fb9ea": "pd-1.0",
    "985c01ca-d9d2-4df5-a8b9-8a6da251d7d4": "migrate-3.7.2",
    "f167eac0-2a23-4e74-bb1c-abdfb5f74a92": "import-sequences",
    "4cfcd389-ed9e-4174-aa99-150f176e8eec": "import-msa",
    "caaca280-2290-4625-b5c0-76bcfb06e9ac": "import-phylotree"
}
# "15aa399f-dd58-433f-8e94-5b2222cd06c9"
# "5b7e9e40-040b-40fc-9db3-7d707fe9617f"
# "8fac3ce8-8796-445f-ac27-4baedadeff3b"
# "21879d8f-1c0e-4f71-92a9-88bc6a3aa14b"
# "83077626-cf8c-48d3-854b-a355afdb7df9"
# "fc1fb247-6b76-420c-9c48-f69f154cbe1d"
# "91a5b4a7-3359-4eed-98df-497c42d0c3c1"
# "6ba7c06b-c164-4049-8259-f713920284a2"
# "7d26e0fe-501e-4568-92a9-6c250eceb705"
# "98da9069-5c62-44c3-8f69-985e439d106d"
# "a2aef599-a34d-4290-bde5-14899b70eff1"
# "1a17a8c1-2755-4bfb-b42e-c65aa53800e6"
# "3e7dab12-dd92-49cc-83eb-d4f1c0e4a62d"
# "86a940ec-3327-4c62-8ca3-f4fbf47a62ce"
# "0b827be9-99c4-42d5-ad9a-8e7f3a4ebc39"
# "dc0ae8e7-b05c-4630-95d5-fd823ca3091a"
# "1882a009-a285-45a4-bcd3-bb8e23bab2b6"
# "7b5e6398-3a35-400c-92cc-551687058cd0"
# "ab579882-0060-45ed-b747-5a43b39c8d25"
# "f74ea2b7-8c19-4ad2-b557-138d5dd4a048"
# "c8df0c20-9cd5-499b-92d4-5fb35b5a369a"
# "16159c67-9325-4f0d-b0c5-2f01588612ea"
# "fcaa3b92-7ba8-4068-8e26-9321a341b53f"
# "b25db3a5-9bb0-4838-b9cc-98f0176876af"
# "9957c469-5a6b-496a-b7c5-184d7eb84688"
# "35b66586-60b2-4764-83fc-d4af276b0636"
# "a80cd569-584d-4f3f-8721-dcb25032476b"
# "fa15b3d0-3cb3-4d34-8005-c90fbde78ac6"
# "49e206a0-e5ac-4848-b47d-8340a728a2e2"
# "c8b62e30-7ec5-462a-903e-83ea1485e995"
# "1f648f74-4b91-4a8d-84d1-c2e3d5c3b5a4"
# "671eb4fa-eeca-4ba6-a88e-ddc8d2ec2b47"
# "b6b76a00-6da8-4c65-99c1-620b806c9dbf"
# "41fa0d02-cc02-4562-b657-994172d1919d"
# "f505261b-d6a1-4024-aa98-7819b0f701e5"
# "54d5633e-9f84-4feb-84b5-ec99de74b668
# "5cd9b251-68de-4aec-9677-faf9724113a8
# 5e5a1bfa-85e8-4ff8-9b1d-8180dc908933
# 75f3373b-ee0e-4e13-85da-2045010f3939
# 6bb4cbbe-d6cd-4c2e-bb59-782e3c9e9f6c
# 44363784-d304-4e7e-a507-8bae48598e50
# 8b62f4aa-d32a-4841-89f5-9ed50da44121
# 5a01d289-8534-40c0-9f56-dc116b609afd
# c87f58b6-cb06-4d39-a0b3-72c2705c5ae1
# c55280d0-f916-4401-a1a4-bb26d8179fd7
# 3e0240e8-b978-48a2-8fdd-9f31f4264064
# ce018826-7b20-4b70-b9b3-168c0ba46eec
# 5b315dc5-ad12-4214-bb6a-bf013f0e4b8c
# f4ab2464-0e36-4d11-8189-a017cab360bc
# 8713fd84-d162-49d3-8549-de0393771836
# 4f7d285a-2d41-4115-8fe9-351c3c703910
# 682b851c-d4c1-467d-b5ca-fbbd1f6b7279
# be5d0563-061d-4588-82a4-8dc0b85cabef
# 6148d4bf-5f55-4a8c-ae27-191c3826fb6e
# e9eda02e-fcec-4773-bc62-51071d6f8e4c
# c98751a2-90de-46cb-b11e-e0097ce77f53
# 598114ac-7d2a-4f5e-9486-2c6ef9cbebfc
# 559443bc-328e-4d6e-be78-cb26576df3cc
# f3ed86f6-650a-4bf1-8d3a-ad984c7ac00b
# 753c138a-81a3-44fb-912d-e42fb610a264
# 43ed5dbc-405a-4245-b390-e5e69e04614c
# d57770a2-69fa-4965-8339-7e91dbf14e5e
# bb862933-984d-4032-8d69-86385cbbfabe
# 06f36c7b-8195-4635-8ed1-a1d63ea864ca
# 4ea0dead-eb5d-4fc9-aa60-fc5d1695630a
# 9bca9dac-a939-4449-ad60-cc9ad995c89a
# 5a3d7582-c4d1-4929-89c7-aa6ffdffe691
# 6d0f63a9-847c-4fdf-93e6-13f8129b6fe8
# 2e461eec-532d-416e-8e93-6aadb897057b
# 0baedab0-5cbf-4dd6-bda1-72bae4f6e067
# cb8c6422-3917-4b0e-8eff-95bca237498e
# 5502986b-d22a-4336-88f5-1397c7ca28f4
# 5ca2b427-39a5-4808-87cc-8a53a540da29
# 666ae4ee-93bc-4257-8cad-3aaffacc628b
# 3e5e94d1-dfcd-46ae-b6b3-b8ed51892bc0
# 003bef41-0a73-4be6-a974-a618e005805b
# e607db06-d036-429c-a6be-7b1c71a07c3e
# d89f8fa6-d07d-4018-91ef-6516c87edf5c
# 672b737e-0247-4d98-b0e7-746f60e7372a
# 5f02f994-87f4-4eeb-924c-8cf79ab00284
# 78f2fe8f-db7d-4894-8c17-e03fc96c36b6
# 8b5a8990-aae0-4090-89bf-2274b8332cb6
# 7b09dd78-f8cd-4942-82d4-c5482b3d417d
# c82f14ef-7fc2-45d0-baf6-034dfc33ebc2
# 47728a0d-f7e4-49be-9fe7-d675f8981514
# b0da1326-a697-4a47-a959-2e779c7c179f
# 087341a7-a371-4b7b-b3b9-92f91825e880
# 25932546-d26c-4367-8c81-0c682094d117
# ec40143f-ae32-4dac-9cfb-caa047e1adb1
# 28615331-4c80-4387-b353-a6fd0338b475
# f7c2c088-5ca1-46a7-90b3-9f446b706724


def initialize_database_data():
    # Load base tables
    load_table(DBSession, Identity, tm_default_users)
    load_table(DBSession, Authenticator, tm_authenticators)
    load_table(DBSession, ObjectType, tm_object_types)
    load_table(DBSession, PermissionType, tm_permissions)
    load_table(DBSession, JobStatus, tm_job_statuses)
    load_table(DBSession, JobManagementType, tm_job_mgmt_types)
    load_table(DBSession, Process, tm_processes)
    # Load a ComputeResource if it does not exist
    cr_id = "f0952819-7a7e-4ed0-a33b-139267fe33f4"
    session = DBSession()
    i = session.query(ComputeResource).filter(ComputeResource.uuid == cr_id).first()
    if not i:
        cr = ComputeResource()
        cr.uuid = cr_id
        cr.name = "balder - ITC"
        jm = session.query(JobManagementType).filter(JobManagementType.name == "galaxy").first()
        cr.jm_type = jm
        cr.jm_location = {"url": "http://balder"}
        cr.jm_credentials = {"user": "", "password": ""}
        session.add(cr)
        session.commit()
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


def register_api(bp: Blueprint, view, endpoint: str, url: str, pk='id', pk_type='int'):
    view_func = view.as_view(endpoint)
    bp.add_url_rule(url, defaults={pk: None}, view_func=view_func, methods=['GET'])
    bp.add_url_rule(url, view_func=view_func, methods=['POST'])
    bp.add_url_rule(f'{url}<{pk_type}:{pk}>', view_func=view_func, methods=['GET', 'PUT', 'DELETE'])
