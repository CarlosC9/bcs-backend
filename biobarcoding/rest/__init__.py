import logging
import os
import sys
from enum import Enum

import redis
import sqlalchemy
from alchemyjsonschema import SchemaFactory, StructuralWalker
from attr import attrs, attrib
from flask import Response, Blueprint, g, request
from flask.views import MethodView
from sqlalchemy import orm, and_, or_
from sqlalchemy.pool import StaticPool
from typing import Dict

import biobarcoding
from biobarcoding.authentication import bcs_session
from biobarcoding.common import generate_json
from biobarcoding.common.helpers import get_module_logger
from biobarcoding.common.pg_helpers import create_pg_database_engine, load_table, load_many_to_many_table, \
    load_computing_resources, load_processes_in_computing_resources, load_process_input_schema, load_table_extended
from biobarcoding.db_models import DBSession, ORMBase, DBSessionChado, ORMBaseChado, ObjectType
from biobarcoding.db_models.bioinformatics import *
from biobarcoding.db_models.sysadmin import *
from biobarcoding.db_models.geographics import *
from biobarcoding.db_models.metadata import *
from biobarcoding.db_models.jobs import *
from biobarcoding.db_models.hierarchies import *

bcs_api_base = "/api"  # Base for all RESTful calls
bcs_gui_base = "/gui"  # Base for the Angular2 GUI
bcs_external_gui_base = "/gui_external"  # Base for the Angular2 GUI when loaded from URL

logger = get_module_logger(__name__)
log_level = logging.DEBUG


def build_json_response(obj, status=200):
    return Response(generate_json(obj),
                    mimetype="text/json",
                    status=status)


class IType(Enum):
    INFO = 1
    WARNING = 2
    ERROR = 3


@attrs
class IssueLocation:
    # TODO IssueLocation
    def __str__(self):
        return f"Put your Ad here :)"


@attrs
class Issue:
    itype = attrib()
    message = attrib()
    location = attrib(default=None)

    def as_dict(self):
        return {"type": self.itype.name, "message": self.message, "location": self.location}

    def __str__(self):
        return f'(type="{self.itype.name}", message="{self.message}", location="{self.location}")'


@attrs
class ResponseObject:
    content = attrib(default=None)  # type: object
    count = attrib(default=0)  # type: int
    issues = attrib(default=[])  # type: List[Issue]
    # Mimetype.
    content_type = attrib(default="text/json")  # type: str
    # HTTTP response status code
    status = attrib(default=200)  # type: int

    def get_response(self) -> Response:
        """
        status == 200. ctype not JSON -> content, ctype, 200
                       ctype     JSON -> dict(issues=issues, content=content), ctype, 200
        status != 200. ctype not JSON -> dict(issues=issues), "text/json", status
                       ctype     JSON -> dict(issues=issues)
        :return:
        """
        if self.status == 200:
            if self.content_type in ("text/json", "application/json"):
                obj = generate_json(dict(issues=[s.as_dict() for s in self.issues],
                                         content=self.content, count=self.count))
            else:
                obj = self.content
        else:
            self.content_type = "application/json"
            obj = generate_json(dict(issues=[s.as_dict() for s in self.issues]))

        return Response(obj, mimetype=self.content_type, status=self.status)


def prepare_default_configuration(create_directories):
    def default_directories(path, tmp_path):
        REDIS_HOST = "redis"
        REDIS_PORT = 6379
        BROKER_URL = f"redis://{REDIS_HOST}:{REDIS_PORT}/0"
        BACKEND_URL = BROKER_URL
        return {
            "DB_CONNECTION_STRING": f'postgres://postgres:postgres@localhost:5432/',
            "CACHE_FILE_LOCATION": f"{tmp_path}/cache",
            "CELERY_BROKER_URL": BROKER_URL,
            "CELERY_BACKEND_URL": BACKEND_URL,
            "REDIS_HOST_FILESYSTEM_DIR": f"{tmp_path}/sessions",
            "GOOGLE_APPLICATION_CREDENTIALS": f"{path}/firebase-key.json",
            "CHADO_DATABASE": "postgres",
            "CHADO_HOST": "localhost",
            "CHADO_PORT": "5432",
            "CHADO_USER": "postgres",
            "CHADO_PASSWORD": "postgres",
            "CHADO_SCHEMA": "public",
            "GALAXY_API_KEY": "fakekey",
            "GALAXY_LOCATION": "http://localhost:8080"
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
    try:
        _, file_name = prepare_default_configuration(False)
        if not os.environ.get(biobarcoding.config_file_var):
            found = False
            logger.debug(f"Trying {file_name}")
            for f in [file_name]:
                if os.path.isfile(f):
                    print(f"Assuming {f} as configuration file")
                    logger.debug(f"Assuming {f} as configuration file")
                    found = True
                    os.environ[biobarcoding.config_file_var] = f
                    break
                else:
                    logger.debug(f"Creating {f} as configuration file")
            if not found:
                cfg, file_name = prepare_default_configuration(True)
                print(f"Generating {file_name} as configuration file:\n{cfg}")
                with open(file_name, "wt") as f:
                    f.write(cfg)
                os.environ[biobarcoding.config_file_var] = file_name

        print("-----------------------------------------------")
        print(f'Configuration file at: {os.environ[biobarcoding.config_file_var]}')
        print("-----------------------------------------------")
        flask_app.config.from_envvar(biobarcoding.config_file_var)
        logger.debug(flask_app.config)

    except Exception as e:
        print(f"{biobarcoding.config_file_var} environment variable not defined!")
        print(e)


def  construct_session_persistence_backend(flask_app):
    # A REDIS instance needs to be available. Check it
    # A local REDIS could be as simple as:
    #
    # docker run --rm -p 6379:6379 redis:alpine
    #
    d = {}
    if 'REDIS_HOST' in flask_app.config:
        r_host = flask_app.config['REDIS_HOST']
        d["SESSION_KEY_PREFIX"] = "bcs:"
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
                print("Trying connection to REDIS '" + r_host + "'")
                rs2.ping()
                print("Connected to REDIS instance '" + r_host + "'")
            except:
                print("REDIS instance '" + r_host + "' not reachable, exiting now!")
                sys.exit(1)
        elif "SESSION_TYPE" not in d:
            print("No session persistence backend configured, exiting now!")
            sys.exit(1)
    return d


# Initialization of tables

tm_object_type_fields = ["id", "uuid", "name"]
tm_object_types = [  # ObjectType
    (bio_object_type_id["sequence"], "22bc2577-883a-4408-bba9-fcace20c0fc8", "sequence"),
    (bio_object_type_id["multiple-sequence-alignment"], "e80a7d27-3ec8-4aa1-b49c-5498e0f85bee",
     "multiple-sequence-alignment"),
    (bio_object_type_id["phylogenetic-tree"], "d30120f0-28df-4bca-90e4-4f0676d1c874", "phylogenetic-tree"),
    (bio_object_type_id["geographic-layer"], "83084df6-7ad0-45d7-b3f1-6de594c78611", "geographic-layer"),
    (bio_object_type_id["sequence-similarity"], "7e23991b-24a0-4da1-8251-c3c3434dfb87", "sequence-similarity"),
    (bio_object_type_id["dataframe"], "5e5a1bfa-85e8-4ff8-9b1d-8180dc908933", "dataframe"),
    (1000, "8fac3ce8-8796-445f-ac27-4baedadeff3b", "sys-functions"),
    # Not bio algorithms, but BCS functions (both backend and frontend), like "export sequence", "browse alignments", "launch bioalgorithm", etc.
    (1001, "21879d8f-1c0e-4f71-92a9-88bc6a3aa14b", "compute-resource"),
    (1002, "83077626-cf8c-48d3-854b-a355afdb7df9", "process")  # Bioinformatic algorithms
]

tm_permissions = {  # PermissionType
    "f19ad19f-0a74-44e8-bd4e-4762404a35aa": "read",
    "04cac7ca-a90b-4d12-a966-d8b0c77fca70": "annotate",
    "d0924822-32fa-4456-8143-0fd48da33fd7": "contribute",
    "83d837ab-01b2-4260-821b-8c4a3c52e9ab": "share",
    "91a5b4a7-3359-4eed-98df-497c42d0c3c1": "execute",
    "d3137471-84a0-4bcf-8dd8-16387ea46a30": "delete"
}

tm_default_users = {  # Identities
    "f3848599-4aa3-4964-b7e1-415d478560be": "admin",
    "2a512320-cef7-41c6-a141-8380d900761b": "_anonymous",
    "27c6a285-dd80-44d3-9493-3e390092d301": "test_user",
    "75f3373b-ee0e-4e13-85da-2045010f3939": "celery_user"
}

tm_default_groups = {
    "f74ea2b7-8c19-4ad2-b557-138d5dd4a048": "all-identified",  # All except anonymous
    "1a17a8c1-2755-4bfb-b42e-c65aa53800e6": "sequence-lab",
    "3e7dab12-dd92-49cc-83eb-d4f1c0e4a62d": "endemisms",
    "86a940ec-3327-4c62-8ca3-f4fbf47a62ce": "land-planning"
}

tm_default_roles = {
    "6ba7c06b-c164-4049-8259-f713920284a2": "sysadmin",
    "7d26e0fe-501e-4568-92a9-6c250eceb705": "data-curator",
    "98da9069-5c62-44c3-8f69-985e439d106d": "data-importer",
    "a2aef599-a34d-4290-bde5-14899b70eff1": "bioinformatic-job-executor",
    "0b827be9-99c4-42d5-ad9a-8e7f3a4ebc39": "basic-user"
    # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx": "celery"
}

tm_authenticators = {  # Authenticator
    "5b7e9e40-040b-40fc-9db3-7d707fe9617f": "firebase",
    "b33193c3-63b9-49f7-b888-ceba619d2812": "firebase-google",
    "c09fa36b-62a3-4904-9600-e3bb5028d809": "firebase-facebook",
    "f510cb30-7a44-4cb1-86f5-1b112e43293a": "firebase-mail",
    "5f32a593-306f-4b69-983c-0a5680556fae": "local",
    "15aa399f-dd58-433f-8e94-5b2222cd06c9": "local-api-key"
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
    "0292821a-dd33-450a-bdd8-813b2b95c456": "ssh",
    "fc1fb247-6b76-420c-9c48-f69f154cbe1d": "ebi"
}

galaxy_tm_processes = {# Preloaded Galaxy processes
    "02f44e54-f139-4ea0-a1bf-fe27054c0d6c": "klustal-1",
    "903a73a9-5a4e-4cec-b8fa-4fc9bd5ffab5": "blast",
    "5c4ba6db-e7f2-4d5c-a89a-76059ac116b1": "mrbayes",
    "ea647c4e-2063-4246-bd9a-42f6a57fb9ea": "pd-1.0",
    "985c01ca-d9d2-4df5-a8b9-8a6da251d7d4": "migrate-3.7.2",
    "f167eac0-2a23-4e74-bb1c-abdfb5f74a92": "import-sequences",
    "4cfcd389-ed9e-4174-aa99-150f176e8eec": "import-msa",
    "caaca280-2290-4625-b5c0-76bcfb06e9ac": "import-phylotree",
    "15aa399f-dd58-433f-8e94-5b2222cd06c9": "Clustal Omega",
    "c8df0c20-9cd5-499b-92d4-5fb35b5a369a": "MSA ClustalW",
    "ec40143f-ae32-4dac-9cfb-caa047e1adb1": "ClustalW-PhyMl",
}

ssh_tm_processes = {# Preloaded SSH processes
    "25932546-d26c-4367-8c81-0c682094d117": "SSHTestProcess"
}

tm_processes = {**ssh_tm_processes, **galaxy_tm_processes}

tm_system_functions = {
    "6ba7c06b-c164-4049-8259-f713920284a2": "get sequence",
    "7d26e0fe-501e-4568-92a9-6c250eceb705": "post sequence",
    "98da9069-5c62-44c3-8f69-985e439d106d": "put sequence",
    "a2aef599-a34d-4290-bde5-14899b70eff1": "delete sequence",
    "1a17a8c1-2755-4bfb-b42e-c65aa53800e6": "get sequences",
    "3e7dab12-dd92-49cc-83eb-d4f1c0e4a62d": "post sequences",
    "86a940ec-3327-4c62-8ca3-f4fbf47a62ce": "put sequences",
    "0b827be9-99c4-42d5-ad9a-8e7f3a4ebc39": "delete sequences",
    "dc0ae8e7-b05c-4630-95d5-fd823ca3091a": "get jobs",
    "1882a009-a285-45a4-bcd3-bb8e23bab2b6": "get job",
    "7b5e6398-3a35-400c-92cc-551687058cd0": "post job",
    "ab579882-0060-45ed-b747-5a43b39c8d25": "delete job"
    # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx": "put job"

}


#
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

def initialize_database_data():
    # Load base tables
    load_table(DBSession, Identity, tm_default_users)
    load_table(DBSession, Authenticator, tm_authenticators)
    load_table_extended(DBSession, ObjectType, tm_object_type_fields, tm_object_types)
    load_table(DBSession, SystemFunction, tm_system_functions)
    load_table(DBSession, PermissionType, tm_permissions)
    load_table(DBSession, JobStatus, tm_job_statuses)
    load_table(DBSession, JobManagementType, tm_job_mgmt_types)
    load_table(DBSession, Process, tm_processes)
    load_table(DBSession, Group, tm_default_groups)
    load_table(DBSession, Role, tm_default_roles)
    load_computing_resources(DBSession)
    load_processes_in_computing_resources(DBSession)
    load_process_input_schema(DBSession)

    # Load default authentication for "test_user"
    session = DBSession()
    test_user_id = "test_user"
    authenticator_id = "5f32a593-306f-4b69-983c-0a5680556fae"  # "local". Just the user name is good to be authorized
    iden = session.query(Identity).filter(Identity.name == test_user_id).first()
    authentication = session.query(Authenticator).filter(Authenticator.uuid == authenticator_id).first()
    iden_authentication = session.query(IdentityAuthenticator). \
        filter(
        and_(IdentityAuthenticator.identity == iden, IdentityAuthenticator.authenticator == authentication)).first()
    if not iden_authentication:
        iden_authentication = IdentityAuthenticator()
        iden_authentication.identity = iden
        iden_authentication.authenticator = authentication
        iden_authentication.name = "test_user"
        iden_authentication.email = "test@test.org"
        session.add(iden_authentication)
    session.commit()

    # Set test_user roles and groups
    load_many_to_many_table(DBSession, RoleIdentity, Role, Identity, ["role_id", "identity_id"],
                            [["sysadmin", "test_user"]])
    load_many_to_many_table(DBSession, GroupIdentity, Group, Identity, ["group_id", "identity_id"],
                            [["all-identified", "test_user"]])

    # Insert a test BOS
    # session = DBSession()
    # seq = Sequence()
    # seq.name = "test seq"
    # session.add(seq)
    # msa = MultipleSequenceAlignment()
    # msa.name = "test msa"
    # session.add(msa)
    # phyt = PhylogeneticTree()
    # phyt.name = "test phylo tree"
    # session.add(phyt)
    # seq_sim = SequenceSimilarity()
    # seq_sim.name = "test seq simil"
    # session.add(seq_sim)
    # geo_layer = GeographicLayer()
    # geo_layer.name = "test geo layer"
    # session.add(geo_layer)
    # session.commit()
    # DBSession.remove()


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
    cfg = flask_app.config
    if 'CHADO_CONNECTION_STRING' in flask_app.config:
        db_connection_string = cfg["CHADO_CONNECTION_STRING"]
    elif 'CHADO_DATABASE' in flask_app.config:
        db_connection_string = f'postgres://{cfg["CHADO_USER"]}:{cfg["CHADO_PASSWORD"]}@{cfg["CHADO_HOST"]}:{cfg["CHADO_PORT"]}/{cfg["CHADO_DATABASE"]}'
        print("Connecting to Chado database server")
        print(db_connection_string)
        print("-----------------------------")
        biobarcoding.chado_engine = sqlalchemy.create_engine(db_connection_string, echo=True)
        # global DBSessionChado # global DBSessionChado registry to get the scoped_session
        DBSessionChado.configure(
            bind=biobarcoding.chado_engine)  # reconfigure the sessionmaker used by this scoped_session
        orm.configure_mappers()  # Important for SQLAlchemy-Continuum
        ORMBaseChado.metadata.bind = biobarcoding.chado_engine
        ORMBaseChado.metadata.reflect()
    else:
        print("No CHADO connection defined (CHADO_CONNECTION_STRING), exiting now!")
        print(flask_app.config)
        sys.exit(1)


from sqlalchemy import text, Table, Column, BigInteger, ForeignKey, Index, MetaData, UniqueConstraint
import subprocess


def initialize_chado_edam(flask_app):
    cfg = flask_app.config
    if 'CHADO_DATABASE' in flask_app.config:
        chado_xml_folder = os.path.join(os.sep, "app", "docker_init", "chado_xml")
        db_relationship_comment = """Specifies relationships between databases.  This is 
        particularly useful for ontologies that use multiple prefix IDs for it\'s vocabularies. For example,
        the EDAM ontology uses the prefixes \"data\", \"format\", \"operation\" and others. Each of these would
        have a record in the db table.  An \"EDAM\" record could be added for the entire ontology to the
        db table and the previous records could be linked as \"part_of\" EDAM.  As another example
        databases housing cross-references may have sub databases such as NCBI (e.g. Taxonomy, SRA, etc).
        This table can use a \'part_of\' record to link all of them to NCBI."""

        db_connection_string = f'postgres://{cfg["CHADO_USER"]}:{cfg["CHADO_PASSWORD"]}@{cfg["CHADO_HOST"]}:{cfg["CHADO_PORT"]}/{cfg["CHADO_DATABASE"]}'
        print("Connecting to Chado database server to insert EDAM")
        print(db_connection_string)
        print("-----------------------------")
        biobarcoding.chado_engine = sqlalchemy.create_engine(db_connection_string, echo=True)

        with biobarcoding.chado_engine.connect() as conn:
            relationship_id = conn.execute(text("select * from INFORMATION_SCHEMA.TABLES where TABLE_NAME = 'db_relationship'")).fetchone()
            if relationship_id is None:
                try:
                    meta = MetaData()
                    meta.reflect(biobarcoding.chado_engine, only=["db",
                                                                  "cvterm"])  # this get information about the db and cvterm tables in the existing database
                    # CREATE db_relationship table to relate EDAM submodules with EDAM itself
                    db_relationship = Table('db_relationship', meta,
                                            Column('db_relationship_id', BigInteger, primary_key=True),
                                            Column('type_id',
                                                   BigInteger,
                                                   ForeignKey('cvterm.cvterm_id', ondelete="CASCADE", initially="DEFERRED"),
                                                   nullable=False),
                                            Column('subject_id',
                                                   BigInteger,
                                                   ForeignKey('db.db_id', ondelete="CASCADE", initially="DEFERRED"),
                                                   nullable=False),
                                            Column('object_id',
                                                   BigInteger,
                                                   ForeignKey('db.db_id', ondelete="CASCADE", initially="DEFERRED"),
                                                   nullable=False),
                                            UniqueConstraint('type_id', 'subject_id', 'object_id',
                                                             name='db_relationship_c1'),
                                            Index('db_relationship_idx1', 'type_id'),
                                            Index('db_relationship_idx2', 'subject_id'),
                                            Index('db_relationship_idx3', 'object_id'),
                                            comment=db_relationship_comment
                                            )

                    # DELETE CVTERMS THAT ARE IN THE CHADO XML AND ALSO IN THE DATABASE
                    conn.execute(text("delete from cvterm where name = 'comment'and cv_id = 6 and is_obsolete = 0"))
                    conn.execute(text("delete from cvterm where name = 'has_function'and cv_id = 4 and is_obsolete = 0"))
                    conn.execute(text("delete from cvterm where name = 'has_input'and cv_id = 4 and is_obsolete = 0"))
                    conn.execute(text("delete from cvterm where name = 'has_output'and cv_id = 4 and is_obsolete = 0"))
                    conn.execute(text("delete from cvterm where name = 'is_a'and cv_id = 4 and is_obsolete = 0"))
                    conn.execute(text("delete from cvterm where name = 'is_anonymous'and cv_id = 6 and is_obsolete = 0"))

                    db_relationship.create(bind=conn)

                except Exception as e:
                    print("Something went wrong when inserting EDAM:")
                    print(e)
            else:
                print("The preprocess to insert EDAM was already executed")

            try:
                edam_id = conn.execute(text("select db_id from db where name=\'EDAM\'")).fetchone()
                if edam_id is None:
                    # INSERT THE CHADO XML (EDAM) IN THE DATABASE
                    command = ["stag-storenode.pl", "-d",
                               f'dbi:Pg:dbname={cfg["CHADO_DATABASE"]};host={cfg["CHADO_HOST"]};port={cfg["CHADO_PORT"]}',
                               "--user", cfg["CHADO_USER"], "--password", cfg["CHADO_PASSWORD"],
                               os.path.join(chado_xml_folder, "EDAM.chado")]

                    subprocess.run(command)
                    #We commit 2 times because we need a commit on the last block to execute this one
                    with biobarcoding.chado_engine.connect() as conn:
                        # FILL THE db_relationship TABLE
                        edam_id = conn.execute(text("select db_id from db where name=\'EDAM\'")).fetchone()[0]
                        db_ids = conn.execute(text("select db_id from db where db_id > " + str(edam_id))).fetchall()
                        type_id = conn.execute(text(
                            "select cvterm_id from cvterm join cv on cvterm.cv_id = cv.cv_id where cvterm.name = \'part_of\' "
                            "and cv.name = \'relationship\'")).fetchone()[0]

                        for db_id in db_ids:
                            conn.execute(text("INSERT INTO db_relationship (type_id,subject_id,object_id)" +
                                              "VALUES(" + str(type_id) + "," + str(db_id[0]) + "," + str(edam_id) + ")"))
                else:
                    print("EDAM was already inserted")

            except Exception as e:
                print("""WARNING! Something went wrong when inserting EDAM. 
                        You will need to delete the database if you do not want to insert EDAM:""")
                print(e)

    else:
        print("No CHADO connection defined (CHADO_CONNECTION_STRING), exiting now!")
        print(flask_app.config)
        sys.exit(1)


def register_api(bp: Blueprint, view, endpoint: str, url: str, pk='id', pk_type='int'):
    view_func = view.as_view(endpoint)
    bp.add_url_rule(url, defaults={pk: None}, view_func=view_func, methods=['GET'])
    bp.add_url_rule(url, view_func=view_func, methods=['POST'])
    bp.add_url_rule(f'{url}<{pk_type}:{pk}>', view_func=view_func, methods=['GET', 'PUT', 'DELETE'])


def make_simple_rest_crud(entity, entity_name: str, execution_rules: Dict[str, str] = {}):
    """
    Create a CRUD RESTful endpoint

    :param entity: the entity class to manage
    :param entity_name: the name of the entity, in plural
    :param execution_rules: a dictionary of execution rules (according to the decorator "bcs_session") by CRUD method
    :return: the blueprint (to register it) and the class (if needed)
    """

    class CrudAPI(MethodView):

        page: int = None
        page_size: int = None

        @bcs_session(read_only=True, authr=execution_rules.get("r"))
        def get(self, _id=None):  # List or Read
            db = g.bcs_session.db_session
            r = ResponseObject()
            if _id is None:
                # List of all
                query = db.query(entity)
                # TODO Pagination
                self.__check_data(request.args)
                if self.page and self.page_size:
                    query = query.offset((self.page - 1) * self.page_size).limit(self.page_size)
                # TODO Detail of fields
                r.content = query.all()
            else:
                # Detail
                # TODO Detail of fields
                r.content = db.query(entity).filter(entity.id == _id).first()

            return r.get_response()

        @bcs_session(authr=execution_rules.get("c"))
        def post(self):  # Create
            db = g.bcs_session.db_session
            r = ResponseObject()
            t = request.json
            s = entity.Schema().load(t, instance=entity())
            db.add(s)
            return r.get_response()

        @bcs_session(authr=execution_rules.get("u"))
        def put(self, _id):  # Update (total or partial)
            db = g.bcs_session.db_session
            r = ResponseObject()
            t = request.json
            s = db.query(entity).filter(entity.id == _id).first()
            s = entity.Schema().load(t, instance=s)
            db.add(s)
            return r.get_response()

        @bcs_session(authr=execution_rules.get("d"))
        def delete(self, _id):  # Delete
            db = g.bcs_session.db_session
            r = ResponseObject()
            s = db.query(Identity).filter(Identity.id == _id).first()
            db.delete(s)
            return r.get_response()

        def __check_data(self, data):
            self.page = int(data.get(self.page, 1))
            self.page_size = int(data.get(self.page_size, 1000000))

    # If the following line is uncommented, it complains on "overwriting an existing endpoint function". This function is public, because it is just the schema, so, no problem.
    # @bcs_session()
    def get_entities_json_schema():
        r = ResponseObject()
        factory = SchemaFactory(StructuralWalker)
        r.content = factory(entity)

        return r.get_response()

    @bcs_session()
    def get_entity_json_schema(_id):
        db = g.bcs_session.db_session
        r = ResponseObject()
        factory = SchemaFactory(StructuralWalker)
        ent = db.query(entity).filter(entity.id == _id).first()
        r.content = dict(schema=factory(entity), data=ent)
        return r.get_response()

    bp_entity = Blueprint(f'bp_{entity_name}', __name__)
    register_api(bp_entity, CrudAPI, entity_name, f"{bcs_api_base}/{entity_name}/", "_id")
    bp_entity.add_url_rule(f"{bcs_api_base}/{entity_name}.schema.json", view_func=get_entities_json_schema,
                           methods=['GET'])
    bp_entity.add_url_rule(f"{bcs_api_base}/{entity_name}/<int:_id>.schema.json", view_func=get_entity_json_schema,
                           methods=['GET'])

    return bp_entity, CrudAPI


# def filter_parse(filter_str: str) -> filter_chado, filter_bcs:
#     """
#     [{"campo1": "<condicion>", "campo2": "<condicion>"}, {"campo1": "<condicion"}]
#     <condicion>: {"op": "<operador", "left": "<valor>", "right": <valor>, "unary": <valor>}
#
#     :param filter_str:
#     :return:
#     """
#     def append_bcs_condition(bcs_and_clause, field, condition):
#         if field == "...":
#             obj = ...
#         elif field == "...":
#             obj = ...
#         op = condition["op"]
#         if condition == "in":
#             v = condition["unary"]
#             cond = obj.in_(v)
#         elif condition == "eq":
#             cond = obj == v
#         elif condition == "between":
#             left = condition["left"]
#             right = condition["right"]
#             cond = obj.between_(left, right)
#
#     filter = json.loads(filter_str)
#     bcs_where = ...
#     chado_where = ...
#     for and_clause in filter:
#         bcs_and_clause = ...
#         chado_and_clause = ...
#         for field, condition in and_clause.items():
#             if field in (...): # Chado fields
#                 append_chado_condition(chado_and_clause, field, condition)
#             else: # BCS
#                 append_bcs_condition(bcs_and_clause, field, condition)
#         concatenate_
#     return filter(bcs_where), filter(chado_where)

def filter_parse(orm, filter, aux_filter=None):
    """
    @param orm: <ORMSqlAlchemy>
    @param filter:
        [{"campo1": "<condicion>", "campo2": "<condicion>"}, {"campo1": "<condicion"}]
        <condicion>: {"op": "<operador", "left": "<valor>", "right": <valor>, "unary": <valor>}
    @param aux_filter: <callable function(filter)>
    @return: <orm_clause_filter>
    """
    def get_condition(orm, field, condition):
        obj = getattr(orm, field)
        op = condition["op"]
        value, left, right = condition.get("unary"), condition.get("left"), condition.get("right")
        if op == "in":
            return obj.in_(value)
        elif op == "eq":
            return obj == value
        elif op == "between":
            return obj.between_(left, right)
        return True

    try:
        if not isinstance(filter, (list, tuple)):
            filter=[filter]
        or_clause = []
        for clause in filter:
            and_clause = []
            for field, condition in clause.items():
                if hasattr(orm, field):
                    and_clause.append(get_condition(orm, field, condition))
                else:
                    print(f'Unknown column "{field}" for the given tables.')
            if aux_filter:
                and_clause+=aux_filter(clause)
            or_clause.append(and_(*and_clause))
        return or_(*or_clause)
    except Exception as e:
        print(e)
        return False


def paginator(query, pagination):
    if 'pageIndex' in pagination and 'pageSize' in pagination:
        page = pagination.get('pageIndex')
        page_size = pagination.get('pageSize')
        return query\
            .offset((page - 1) * page_size)\
            .limit(page_size)
    return query