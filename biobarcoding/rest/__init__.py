import json
import logging
import os
import subprocess
import traceback

import sys
from enum import Enum
from typing import Dict, List, Callable
from urllib.parse import unquote

import redis
from alchemyjsonschema import SchemaFactory, StructuralWalker
from attr import attrs, attrib
from flask import Response, Blueprint, g, request
from flask.views import MethodView
import sqlalchemy
from sqlalchemy import orm, and_, or_
from sqlalchemy.pool import StaticPool
from bioblend import galaxy

from .. import config_file_var, app_acronym
from ..authentication import n_session
from ..common import generate_json, ROOT
from ..common.pg_helpers import create_pg_database_engine, load_table, load_many_to_many_table, \
    load_computing_resources, load_processes_in_computing_resources, load_process_input_schema, load_table_extended
from ..db_models import DBSession, DBSessionChado, ORMBaseChado, DBSessionGeo
from ..db_models.bioinformatics import *
from ..db_models.geographics import *
from ..db_models.hierarchies import HierarchyType, Hierarchy, HierarchyNode
from ..db_models.jobs import *
from ..db_models.sa_annotations import *
from ..db_models.sysadmin import *
from ..rest.socket_service import SocketService
import biobarcoding

app_api_base = "/api"  # Base for all RESTful calls
app_gui_base = "/gui"  # Base for the Angular2 GUI
app_external_gui_base = "/gui_external"  # Base for the Angular2 GUI when loaded from URL
app_proxy_base = "/pxy"  # Base for the Backend Proxy

log_level = logging.DEBUG
logger = logging.getLogger('app_logger')

logger.setLevel(log_level)  # set logger level
logFormatter = logging.Formatter("%(name)-12s %(asctime)s %(levelname)-8s %(filename)s:%(funcName)s %(message)s")
consoleHandler = logging.StreamHandler(sys.stdout)  # set streamhandler to stdout
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)


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
    # HTTP response status code
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


def get_default_configuration_dict():
    from appdirs import AppDirs
    app_dirs = AppDirs(f"{app_acronym}-backend")
    # Default directories, multi-platform
    data_path = app_dirs.user_data_dir
    cache_path = app_dirs.user_cache_dir

    REDIS_HOST = "redis"
    REDIS_PORT = 6379
    BROKER_URL = f"redis://{REDIS_HOST}:{REDIS_PORT}/0"
    BACKEND_URL = BROKER_URL

    return dict(
        # SYSTEM DB
        DB_CONNECTION_STRING="postgresql://postgres:postgres@localhost:5432/",
        # CHADO (MOLECULAR DATA DB)
        CHADO_DATABASE="postgres",
        CHADO_HOST="localhost",
        CHADO_PORT="5432",
        CHADO_USER="postgres",
        CHADO_PASSWORD="postgres",
        CHADO_SCHEMA="public",
        # CELERY (TASK EXECUTION)
        CELERY_BROKER_URL=f"{BROKER_URL}",
        CELERY_BACKEND_URL=f"{BACKEND_URL}",
        REDIS_HOST="filesystem:local_session",
        REDIS_HOST_FILESYSTEM_DIR=f"{cache_path}/sessions",
        # GEO (GEOSPATIAL DATA)
        GEOSERVER_USER="admin",
        GEOSERVER_PASSWORD=f"{app_acronym}_ad37",
        GEOSERVER_HOST="geoserver",
        GEOSERVER_PORT="8080",
        # GEOSERVER address from App-Backend, used by the GeoserverProxy
        GEOSERVER_URL="localhost:9180",
        # PostGIS address from App-Backend
        POSTGIS_CONNECTION_STRING="postgresql://postgres:postgres@localhost:5435/",
        # PostGIS address from Geoserver ("host", "port" could differ from those in POSTGIS_CONNECTION_STRING)
        POSTGIS_USER="postgres",
        POSTGIS_PASSWORD="postgres",
        POSTGIS_PORT="5432",
        POSTGIS_HOST="localhost",
        POSTGIS_DB=f"{app_acronym}_geoserver",
        # COMPUTE RESOURCES
        RESOURCES_CONFIG_FILE_PATH=f"{data_path}/compute_resources.json",
        JOBS_LOCAL_WORKSPACE=os.path.expanduser(f'~/{app_acronym}_jobs'),
        SSH_JOBS_DEFAULT_REMOTE_WORKSPACE="/tmp",
        GALAXY_API_KEY="fakekey",
        GALAXY_LOCATION="http://localhost:8080",
        # MISC
        COLOR_PALETTES=f"{data_path}/color_palettes.yaml",
        GOOGLE_APPLICATION_CREDENTIALS=f"{data_path}/firebase-key.json",  # Firebase
        ENDPOINT_URL="http://localhost:5000",  # Self "app-backend" address (so Celery can update things)
        COOKIES_FILE_PATH=f"/tmp/{app_acronym}-cookies.txt",  # Where cookies are stored by Celery
        CACHE_FILE_LOCATION=f"{cache_path}/cache",  # Cached things
        SAMESITE_NONE="True",
        TESTING="True",
        SELF_SCHEMA="",
        ACL_ENABLED="False",
    )


def prepare_default_configuration(create_directories):
    default_cfg = get_default_configuration_dict()

    from appdirs import AppDirs
    app_dirs = AppDirs(f"{app_acronym}-backend")

    # Default directories, multi-platform
    data_path = app_dirs.user_data_dir
    cache_path = app_dirs.user_cache_dir
    if create_directories:
        os.makedirs(data_path, exist_ok=True)
        os.makedirs(cache_path, exist_ok=True)
        directories = [v for k, v in default_cfg.items() if k in ["CACHE_FILE_LOCATION", "REDIS_HOST_FILESYSTEM_DIR"]]
        # Obtain and create directories
        for v in directories:
            os.makedirs(v, exist_ok=True)

    # Default configuration
    return f"""{os.linesep.join([f'{k}="{v}"' for k, v in default_cfg.items()])}""", \
           data_path + os.sep + f"{app_acronym}_local.conf"


def complete_configuration_file(file_name):
    default_cfg = get_default_configuration_dict()
    default_cfg_keys = list(default_cfg.keys())
    with open(file_name, 'r+') as file:
        constants = file.read().splitlines()
        for const in constants:
            if const and not const.strip().startswith('#'):
                key, _ = const.split("=")
                if key.strip() in default_cfg_keys:
                    default_cfg_keys.remove(key.strip())
                else:
                    print(f"WARNING: {key.strip()} is in the configuration file but not in the default config")
        for key in default_cfg_keys:
            print(f"{key} inserted to {file_name} with value {default_cfg[key]}")
            file.write(f'\n{key}="{default_cfg[key]}"')


def load_configuration_file(flask_app):
    # Initialize configuration
    try:
        _, file_name = prepare_default_configuration(False)
        if not os.environ.get(config_file_var):
            logger.debug(f"Trying {file_name}")
            if os.path.isfile(file_name):
                print(f"Assuming {file_name} as configuration file")
                logger.debug(f"Assuming {file_name} as configuration file")
                found = True
                complete_configuration_file(file_name)
                os.environ[config_file_var] = file_name
            else:
                found = False
                logger.debug(f"Creating {file_name} as configuration file")
            if not found:
                cfg, file_name = prepare_default_configuration(True)
                print(f"Generating {file_name} as configuration file:\n{cfg}")
                with open(file_name, "wt") as f:
                    f.write(cfg)
                os.environ[config_file_var] = file_name
        else:
            complete_configuration_file(os.environ.get(config_file_var))
        print("-----------------------------------------------")
        print(f'Configuration file at: {os.environ[config_file_var]}')
        print("-----------------------------------------------")
        flask_app.config.from_envvar(config_file_var)
        logger.debug(flask_app.config)

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
        d["SESSION_KEY_PREFIX"] = f"{app_acronym}:"
        d["SESSION_PERMANENT"] = False
        rs2 = None
        if r_host == "redis_lite":
            try:
                import redislite
                rs2 = redislite.Redis(f"tmp_{app_acronym}_backend_redislite.db")  # serverconfig={'port': '6379'}
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


# SOCKET INITIALIZATION
def init_socket(socketio):
    SocketService(socketio)


# Initialization of tables

tm_object_type_fields = ["id", "uuid", "name"]
tm_object_types = [  # ObjectType
    (data_object_type_id["sequence"], "22bc2577-883a-4408-bba9-fcace20c0fc8", "sequence"),
    (data_object_type_id["multiple-sequence-alignment"], "e80a7d27-3ec8-4aa1-b49c-5498e0f85bee",
     "multiple-sequence-alignment"),
    (data_object_type_id["phylogenetic-tree"], "d30120f0-28df-4bca-90e4-4f0676d1c874", "phylogenetic-tree"),
    (data_object_type_id["sequence-similarity"], "7e23991b-24a0-4da1-8251-c3c3434dfb87", "sequence-similarity"),
    (data_object_type_id["dataframe"], "964f5b72-197a-4046-993b-a54e5896fd24", "dataframe"),
    (data_object_type_id["geolayer"], "b7c404e0-8b5d-4327-b830-85f8209c06e3", "geolayer"),
    (data_object_type_id["grid"], "ffb8e897-5df8-4a01-922c-7901b63ef325", "grid"),
    (data_object_type_id["geoprocess"], "ed678c7d-d4dd-4ddc-b811-a66c030ae574", "geoprocess"),  # CProcess
    (data_object_type_id["geoprocess_instance"], "f0bcbe79-89dc-42fb-b7a8-a103d361f27e", "geoprocess_instance"),
    (data_object_type_id["case_study"], "5f4c1666-e509-436b-8844-082ebe88b2b9", "case_study"),
    (1000, "ad83dcb0-e479-4a44-acf5-387b9731e8da", "sys-functions"),
    (1001, "aa97ddad-8937-4590-afc9-dc00c5601f2a", "compute-resource"),
    (1002, "a4b4a7d2-732f-4db9-9a32-c57c00881eb7", "process")  # Algorithms
]

tm_permissions_fields = ["uuid", "name", "rank"]
tm_permissions = [  # PermissionType
    ("f19ad19f-0a74-44e8-bd4e-4762404a35aa", "read", 1),
    ("04cac7ca-a90b-4d12-a966-d8b0c77fca70", "annotate", 2),
    ("d0924822-32fa-4456-8143-0fd48da33fd7", "contribute", 2),
    ("83d837ab-01b2-4260-821b-8c4a3c52e9ab", "share", 2),
    ("91a5b4a7-3359-4eed-98df-497c42d0c3c1", "execute", 2),
    ("fb1bdaec-3379-4930-b3f8-fea67d0783b7", "modify", 3),
    ("62885891-db3c-4a06-bc28-a407ffdb08a4", "change-permissions", 4),
    ("d3137471-84a0-4bcf-8dd8-16387ea46a30", "delete", 5)
]

tm_default_users = {  # Identities
    "0fc2b361-847c-4b2e-8fbd-533092133eef": "admin",
    "74d20b2c-5b49-462c-8d19-8a72f47b5d1b": "_anonymous",
    "91c8008e-97d5-440c-b3fc-f5a409a44768": "test_user",
    "1c8b0500-32d2-40ce-8da0-4fc772f4c3a3": "celery_user"
}

tm_default_groups = {
    "cad9f7e0-8a91-491c-9012-8900177e132b": "all-identified",  # All except anonymous (not implemented)
    "d62bd35a-b34b-43bb-a359-323d57f95f78": "cc-lab",
    "ff5d8e35-01e8-450c-b0b5-7e4487ba4811": "governance-planning",
    "1a17a8c1-2755-4bfb-b42e-c65aa53800e6": "sequence-lab",
    "3e7dab12-dd92-49cc-83eb-d4f1c0e4a62d": "endemisms",
    "86a940ec-3327-4c62-8ca3-f4fbf47a62ce": "land-planning"
}

tm_default_roles = {
    "c79b4ff7-9576-45f6-a439-551ac23c563b": "sys-admin",
    "6e7ff87f-056b-46b0-b636-c7c0eb48cabe": "owner",
    "0cb4ea40-26a6-4960-9c17-dd3012a53e0f": "basic-user",
    "feb13f20-4223-4602-a195-a3ea14615982": "guest",
    "94f3d8ab-7f20-426d-a881-0c19badd2a3d": "data-curator",
    "6e6b4cce-74f0-44f7-ae82-a799ade98739": "data-importer",
    "356e2030-28f5-49d5-a788-db553a889736": "data-updater",
    "abb2cfe7-7990-4d84-b119-0834cc9937b5": "input-layers-admin",
    "4b9e74b3-b0cc-4861-80ae-8390cd54cfab": "geoprocess-instances-admin",
    "a2aef599-a34d-4290-bde5-14899b70eff1": "bioinformatic-job-executor"
}

tm_authenticators = {  # Authenticator
    "5b7e9e40-040b-40fc-9db3-7d707fe9617f": "firebase",
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
    "65171b3f-e13d-4f5b-8604-2331d292baa5": "cancelling",
    "013b2f3b-4b2f-4b6c-8f5f-425132aea74b": "cancelled",
    "3eef41be-fde1-4ad4-92d0-fe795158b41d": "paused",
}

tm_job_mgmt_types = {
    "38fb34f7-a952-4036-9b0b-4d6c59e8f8d4": "galaxy",
    "0292821a-dd33-450a-bdd8-813b2b95c456": "ssh",
    "fc1fb247-6b76-420c-9c48-f69f154cbe1d": "ebi"
}

tm_processes = {  # Preloaded processes
    "64249aae-52c6-4225-bbcf-be29286d00af": "geoprocess",
    "c8df0c20-9cd5-499b-92d4-5fb35b5a369a": "MSA ClustalW",
    "ec40143f-ae32-4dac-9cfb-caa047e1adb1": "ClustalW-PhyMl",
    "c87f58b6-cb06-4d39-a0b3-72c2705c5ae1": "PAUP Parsimony",
    "3e0240e8-b978-48a2-8fdd-9f31f4264064": "Phylogenetic Diversity Analyzer",
    "c55280d0-f916-4401-a1a4-bb26d8179fd7": "MSA ClustalW + PAUP Parsimony",
    "ce018826-7b20-4b70-b9b3-168c0ba46eec": "PAUP Parsimony + Phylogenetic Diversity Analyzer",
    "5b315dc5-ad12-4214-bb6a-bf013f0e4b8c": "MSA ClustalW + PAUP Parsimony + Phylogenetic Diversity Analyzer",
}

tm_system_functions = {
    # API
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
    "ab579882-0060-45ed-b747-5a43b39c8d25": "delete job",
    # GUI (access to screens)
    "10ad3429-79f4-430e-9474-e757911985f9": "gui-layers",
    "bbd6bb13-7256-4e21-ace5-ce64712e63ce": "gui-layer-import",
    "351b97c7-1b55-480d-bd3f-5090739632a9": "gui-layer-export",
    "ffc1e390-62df-4c1b-b42c-cd1e50cb5e23": "gui-layer-read",
    "b71e08df-1d48-4d64-813a-a476e04e4e68": "gui-layer-edit",
    "2a903da7-0d33-4f41-9652-9fd5ee729f82": "gui-layer-delete",
    "887146f7-4781-45db-8761-592188d03b53": "gui-layer-acl",
    "fd0981d1-e59e-4775-b3f1-8f6a7a70911d": "gui-layers-graphical",
    "040ca696-380a-4859-99b2-1c06b00080f1": "gui-geoprocesses",
    "87bb3e53-ded7-4fcc-b6a1-70422989627f": "gui-geoprocess-read",
    "0d07a7c1-6e91-4948-8b0e-45604b6ed4b4": "gui-geoprocess-edit",
    "2e13521b-5bb1-4d4d-91e7-edf36f0c06bf": "gui-instances-candidate",
    "3b1f8929-218a-403a-bf14-1d33fa5e0782": "gui-instance-read",
    "ede453bb-8203-44a3-a78c-312651036667": "gui-instances-schedule",
    "59e683eb-b4d1-442a-b765-f0f06eaecdc4": "gui-instance-delete",
    "2fe1df59-2d1b-4d43-ae03-f5a25fa33850": "gui-instance-create",
    "570330e9-7748-4b4b-a097-62b22df5aa12": "gui-instance-governance",
    "aab7b4c6-635c-4585-9c52-f672c786854c": "gui-instance-update",
    "6092fab1-382a-4ad7-9e22-fe29bc127b81": "gui-instances-scheduled",
    "45e78495-baa8-4bb3-92b8-02744ba391e6": "gui-instances-finished",
    "dd296392-8dac-43f8-a2c9-e7ae0d37ce0a": "gui-port-types",
    "04031418-c368-4242-bcb3-77fd4ae82f8e": "gui-jobs-executing",
    "52dc01da-4c01-4c32-a6e5-771c6ec93632": "gui-job-read",
    "367ed1f2-56de-4b38-9440-574416b1adbd": "gui-job-cancel",
    "780aaec6-030e-47ff-800a-e5d4b46157c2": "gui-jobs-finished",
    "e15847a3-4b93-47cc-88b1-7eba8982ef4b": "gui-case-studies",
    "b4aa0e27-b88e-44d3-aaa5-f3b62be57d92": "gui-case-study-create",
    "3fe9c93d-9927-4f95-926d-8b891e202b16": "gui-case-study-read",
    "3fb0f0bc-2234-4c61-97d2-45c3a25b61ff": "gui-case-study-edit",
    "e74c05c1-bd45-4c45-b91e-b01fc11468a4": "gui-subjects",
    "24386985-fa3b-4ba6-91ba-333c57082929": "gui-sources",
    "1a8f6e2b-c9c6-414e-a9ee-011f416130db": "gui-crs",
    "e8582e8d-9dd0-46cb-aa07-9ce1d57b2b32": "gui-users",
    "505828db-2028-41d6-a151-212b8084028b": "gui-users-authorize",
    "54cd5ed1-52cf-45b9-b58c-8d256b100dad": "gui-user-read",
    "0816cc73-8fa2-4dc5-bbfc-19826c4899fa": "gui-user-edit",
    "39382b18-6a6b-4a17-8925-ec86048aaa43": "gui-user-delete",
    "b5d05b07-196a-4e93-b1d1-4fec95428307": "gui-roles",
    "440b73f5-fdfa-463b-9bb8-ae4e8c3c3fec": "gui-role-edit",
    "c2d1e08a-3637-432e-8fc0-ccd8990bce6e": "gui-role-delete",
    "533c142e-86de-4a7a-941d-31faab8eb64d": "gui-organizations",
    "8ffa434b-12d9-4689-8b5c-29ce838ab80b": "gui-organization-read",
    "216222b4-1f83-40c6-8cd4-e4faeb543d5e": "gui-organization-edit",
    "b22fd110-73cd-42cb-b771-67fc8cbea1b5": "gui-organization-delete",
    "d7990a07-fe58-43bc-acb9-79f322e86f59": "gui-groups",
    "93fb2a7f-29cc-4846-a02e-615a768a72c9": "gui-authenticators",
    "f57827dd-3e43-4bcf-8a89-0d9256c68b3c": "gui-permission-types",
    "90665142-905c-42af-90dd-40a907403513": "gui-system-functions"
}

tm_system_functions_rules = {
    "gui-layers": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin', 'guest')",
    "gui-layer-import": "role in ('sys-admin', 'input-layers-admin')",
    "gui-layer-export": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin')",
    "gui-layer-read": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin', 'guest')",
    "gui-layer-edit": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin')",
    "gui-layer-delete": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin')",
    "gui-layer-acl": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin')",
    "gui-layers-graphical": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin', 'guest')",
    "gui-geoprocesses": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-geoprocess-read": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-geoprocess-edit": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-instances-candidate": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-instance-read": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-instances-schedule": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-instance-delete": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-instance-create": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-instance-governance": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-instance-update": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-instances-scheduled": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-instances-finished": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-port-types": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-jobs-executing": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-job-read": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-job-cancel": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-jobs-finished": "role in ('sys-admin', 'geoprocess-instances-admin')",
    "gui-case-studies": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin')",
    "gui-case-study-create": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin')",
    "gui-case-study-read": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin')",
    "gui-case-study-edit": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin')",
    "gui-subjects": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin', 'guest')",
    "gui-sources": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin', 'guest')",
    "gui-crs": "role in ('sys-admin', 'input-layers-admin', 'geoprocess-instances-admin', 'guest')",
    "gui-users": "role in ('sys-admin')",
    "gui-users-authorize": "role in ('sys-admin')",
    "gui-user-read": "role in ('sys-admin')",
    "gui-user-edit": "role in ('sys-admin')",
    "gui-user-delete": "role in ('sys-admin')",
    "gui-roles": "role in ('sys-admin')",
    "gui-role-edit": "role in ('sys-admin')",
    "gui-role-delete": "role in ('sys-admin')",
    "gui-organizations": "role in ('sys-admin')",
    "gui-organization-read": "role in ('sys-admin')",
    "gui-organization-edit": "role in ('sys-admin')",
    "gui-organization-delete": "role in ('sys-admin')",
    "gui-groups": "role in ('dont-show')",
    "gui-authenticators": "role in ('dont-show')",
    "gui-permission-types": "role in ('dont-show')",
    "gui-system-functions": "role in ('dont-show')"
}

tm_browser_filter_form_fields = ("id", "uuid")
tm_browser_filter_forms = [
    (data_object_type_id["sequence"], "6bb4cbbe-d6cd-4c2e-bb59-782e3c9e9f6c"),
    (data_object_type_id["multiple-sequence-alignment"], "44363784-d304-4e7e-a507-8bae48598e50"),
    (data_object_type_id["phylogenetic-tree"], "8b62f4aa-d32a-4841-89f5-9ed50da44121"),
    (data_object_type_id["geolayer"], "b4d58b86-114d-47fd-b813-ce3fcf18cbf1"),
    (data_object_type_id["grid"], "dda18cf6-0449-4824-9a6d-3b86760f27ff"),
    (data_object_type_id["dataframe"], "c645a547-49c4-49d9-a151-a81ad404c38e"),
    (100, "a2781f70-46b9-45c5-b3e7-48295d93d339"),  # Processes
    (101, "409e25a3-c2d2-468f-ac29-ab0493df2b5c"),  # Process executions
]

ht_cl_name = "Code List"
tm_hierarchy_types = {
    "45e91fea-6ca3-40c8-a798-98ca6a8e0e2a": ht_cl_name,
    "f3625ee6-99e1-4db9-82e4-adb6c1a01762": "Hierarchical Code List",
}

h_subjects_name = "Temas"
h_sources_name = "Fuentes"
h_crs_name = "CRS"
tm_hierarchies_fields = ((HierarchyType, "name", "id", "h_type_id"), "uuid", "name")
tm_hierarchies = [
    (ht_cl_name, "9975f4a0-a321-409a-bb2b-7afd05a60117", h_subjects_name),
    (ht_cl_name, "8267366d-97c7-47d8-83b5-68fe763c5e81", h_sources_name),
    (ht_cl_name, "7cdc80c7-739d-4787-8529-b52f83994930", h_crs_name)
]

tm_code_list_fields = ((Hierarchy, "name", "id", "hierarchy_id"), "uuid", "name")
tm_code_list_subjects = [
    (h_subjects_name, "2d1d71be-2e7e-4c0a-8d4c-21d2bc5bdad6", "Altimetría y modelos digitales"),
    (h_subjects_name, "eb85c308-61e3-4be6-b2c4-4d86b46a62f8", "Ámbitos de Ordenación Específicos"),
    (h_subjects_name, "68ecc77c-bc78-4352-a407-ee7bb06c04d1", "Áreas de Especial Protección"),
    (h_subjects_name, "08bd7b77-a977-4447-962b-78c0796bc06a", "Cartografía Específica"),
    (h_subjects_name, "dd1cd206-a45c-43d4-95aa-a0f1ac3c5f96", "Catastro"),
    (h_subjects_name, "1247d2cd-dcd1-4238-b44f-ea816c97b084", "Clima"),
    (h_subjects_name, "cb6b0956-4511-42d9-a75b-3b45f094a786", "Coberturas físicas y biológicas"),
    (h_subjects_name, "fdef90ae-c9c1-4415-8bf8-19b31799beb5", "Datos Estadísticos"),
    (h_subjects_name, "343e8354-925f-4c8e-aefa-7834ff2c3cc6", "Edafología"),
    (h_subjects_name, "763a0def-dc78-4561-9f8d-042d09b36ff9", "Elementos hidrográficos"),
    (h_subjects_name, "c5832d73-ba79-46e0-9448-21d0c101927f", "Entidades de población"),
    (h_subjects_name, "efc017c3-45f8-4811-9262-bfb92efc3204", "Equipamientos y servicios"),
    (h_subjects_name, "3fa8c241-05e8-40a3-8f75-e193bfe69999", "Geocodificaciones"),
    (h_subjects_name, "a26d12d3-ead3-4681-95d6-c6fb0c6cbb10", "Infraestructuras y redes de transporte"),
    (h_subjects_name, "1bfd9a77-891d-4863-aa33-09c5976ee823", "Instalaciones agrarias y acuicultura"),
    (h_subjects_name, "57acf475-21d0-4270-a947-5b20ab299539", "Instalaciones Industriales"),
    (h_subjects_name, "f5c2afb6-27ff-4118-91a5-75e6a9541763", "Malla demográfica"),
    (h_subjects_name, "42a78b57-b610-4b66-af4d-ebf090666cba", "Nomenclator de Población"),
    (h_subjects_name, "86ee31c5-bb1a-4c88-af9c-298084920c19", "Ortofotos"),
    (h_subjects_name, "d57653c7-cf75-4333-8677-c58bdf9b0a93", "Otros sin asignación"),
    (h_subjects_name, "809b6d03-db8f-4543-bc4d-93a93dca1839", "Protección Territorial"),
    (h_subjects_name, "8d746d7f-2309-4389-92d1-54ae6015a78f", "Referencias Geográficas"),
    (h_subjects_name, "0afe9d50-7ce8-41ec-9ef2-55712258c43f", "Riesgos"),
    (h_subjects_name, "1497662c-e31a-418c-ae5e-e0eb8e98888e", "Salud y seguridad humana"),
    (h_subjects_name, "ffadcde3-60bf-43e5-9aec-1b54d89370c9", "Unidades Estadísticas"),
    (h_subjects_name, "b8e58143-fe0e-440a-995b-576d2e1e9d51", "Usos del Suelo")
]

tm_code_list_crs = [
    (h_crs_name, "9c534897-2fd6-4165-b161-c0317ff59718", "EPSG:4326"),
    (h_crs_name, "9cbb027f-e532-4d81-bb66-467a06008429", "EPSG:32628"),
]

tm_code_list_sources = [
    (h_sources_name, "d77abee2-17ff-4063-8b54-a8e35647a5f0", "Grafcan"),
    (h_sources_name, "d05c68e7-5674-4b10-8cea-7b9cb8f1b798", "ISTAC"),
    (h_sources_name, "b2be8a7a-67f2-4fd1-afee-cb89164767a9", "ITC"),
    (h_sources_name, "55e9958c-2f73-47c4-afa6-37bf23c683f9", "EUROSTAT"),
]

#
# feb13f20-4223-4602-a195-a3ea14615982
# 5ac8ed2c-20d7-4905-9184-f584947719dd
# 640de355-f58e-4446-a12a-2a99f2cfb2eb
# a31b2a73-d2b2-4d50-a64c-8d6d86559a02
# 7ece7ead-4de9-436f-a12d-86350129e444
# 071a22e8-1d90-4900-936a-0466880480d4
# faec0128-c0e3-4c2c-bbb2-bd852430eed3
# b4f90851-250e-425f-9a94-b677b593c842
# d04fc218-8dea-4938-a4d0-49e3d8174715
# 07c7512b-d7ce-4162-85c4-cda2503f290c
# 95fd8186-927c-4712-8dbc-7fe9dbef35bb
# 848b46b0-8602-42a2-a3fd-1b9be728d729
# 79668ea7-80fa-4327-a433-721a69582542
# 35313fd1-484e-48af-b6e2-e4b7464abb64
# dc50990e-f4ad-4ef3-80c8-f68fd0b1a412
# cd9e5139-8f9b-495f-955d-61577f8531a4
# 611465a5-e57e-4cf7-907f-345b5760d235
# bcaaf1ef-abbf-4b99-8638-5d2bc604998d
# d5a8c860-c6b8-49b3-9447-f75dcb650b1f
# 5b8dd6cc-493d-4994-8fe3-111c83942d0e
# 82dc3955-5125-4341-9fe1-2bc455726bee
# edfb62c9-76dc-4bc7-9282-f927a13711d6
# 94822821-e49b-4011-8188-a7dc026fd094
# 88212405-72cb-49aa-a08f-723ab115b4ec
# 6764cf5b-021a-4922-922a-04754d3c7077
# 6ca59704-dc56-4c91-951e-386802e47dca


def initialize_database_data():
    # Load base tables
    load_table(DBSession, Identity, tm_default_users)
    load_table(DBSession, Authenticator, tm_authenticators)
    load_table_extended(DBSession, ObjectType, tm_object_type_fields, tm_object_types)
    load_table(DBSession, SystemFunction, tm_system_functions)
    load_table_extended(DBSession, PermissionType, tm_permissions_fields, tm_permissions)
    load_table(DBSession, JobStatus, tm_job_statuses)
    load_table(DBSession, JobManagementType, tm_job_mgmt_types)
    load_table(DBSession, Process, tm_processes)
    load_table(DBSession, Group, tm_default_groups)
    load_table(DBSession, Role, tm_default_roles)
    load_computing_resources(DBSession)
    load_processes_in_computing_resources(DBSession)
    load_process_input_schema(DBSession)
    load_table(DBSession, HierarchyType, tm_hierarchy_types)
    load_table_extended(DBSession, Hierarchy, tm_hierarchies_fields, tm_hierarchies)
    load_table_extended(DBSession, HierarchyNode, tm_code_list_fields, tm_code_list_subjects)
    load_table_extended(DBSession, HierarchyNode, tm_code_list_fields, tm_code_list_sources)
    load_table_extended(DBSession, HierarchyNode, tm_code_list_fields, tm_code_list_crs)
    session = DBSession()
    # ACLs for system functions
    system_functions = session.query(SystemFunction).filter(SystemFunction.name.in_(tm_system_functions_rules.keys())).all()
    system_function_uuids = {f.name: f.uuid for f in system_functions}
    obj_type = session.query(ObjectType).filter(ObjectType.name == "sys-functions").first()
    for function, expression in tm_system_functions_rules.items():
        uuid = system_function_uuids[function]
        acl = session.query(ACL).filter(and_(ACL.object_type == obj_type.id, ACL.object_uuid == uuid)).first()
        if not acl:
            acl = ACL()
            acl.object_type = obj_type.id
            acl.object_uuid = uuid
            session.add(acl)
            acl_expression = ACLExpression()
            acl_expression.acl = acl
            acl_expression.expression = expression
            session.add(acl_expression)
        else:
            acl_expression = session.query(ACLExpression).filter(ACLExpression.acl == acl).first()
            if not acl_expression:
                acl_expression = ACLExpression()
                acl_expression.acl = acl
                acl_expression.expression = expression
                session.add(acl_expression)
            else:
                acl_expression.expression = expression
    # load_table_extended(DBSession, BrowserFilterForm, tm_browser_filter_form_fields, tm_browser_filter_forms)

    # Load default authentication for "test_user"
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
        iden_authentication.name = test_user_id
        iden_authentication.email = "test@test.org"
        session.add(iden_authentication)

    # Celery user, also with local authentication
    celery_user_id = "celery_user"
    iden = session.query(Identity).filter(Identity.name == celery_user_id).first()
    iden_authentication = session.query(IdentityAuthenticator). \
        filter(
        and_(IdentityAuthenticator.identity == iden, IdentityAuthenticator.authenticator == authentication)).first()
    if not iden_authentication:
        iden_authentication = IdentityAuthenticator()
        iden_authentication.identity = iden
        iden_authentication.authenticator = authentication
        iden_authentication.name = celery_user_id
        iden_authentication.email = "celery@celery.org"
        session.add(iden_authentication)

    # _anonymous user, also with local authentication
    anonymous_user_id = "_anonymous"
    iden = session.query(Identity).filter(Identity.name == anonymous_user_id).first()
    iden_authentication = session.query(IdentityAuthenticator). \
        filter(
        and_(IdentityAuthenticator.identity == iden, IdentityAuthenticator.authenticator == authentication)).first()
    if not iden_authentication:
        iden_authentication = IdentityAuthenticator()
        iden_authentication.identity = iden
        iden_authentication.authenticator = authentication
        iden_authentication.name = anonymous_user_id
        iden_authentication.email = "_@anonymous.org"
        session.add(iden_authentication)

    session.commit()

    # Set test_user roles and groups
    load_many_to_many_table(DBSession, RoleIdentity, Role, Identity, ["role_id", "identity_id"],
                            [("sys-admin", test_user_id),
                             ("sys-admin", celery_user_id),
                             ("guest", anonymous_user_id)])
    load_many_to_many_table(DBSession, GroupIdentity, Group, Identity, ["group_id", "identity_id"],
                            [("all-identified", test_user_id)])
    load_many_to_many_table(DBSession, ObjectTypePermissionType, ObjectType, PermissionType,
                            ["object_type_id", "permission_type_id"],
                            [["sequence", "read"],
                             ["sequence", "annotate"],
                             ["sequence", "delete"],
                             ["multiple-sequence-alignment", "read"],
                             ["multiple-sequence-alignment", "annotate"],
                             ["multiple-sequence-alignment", "delete"],
                             ["phylogenetic-tree", "read"],
                             ["phylogenetic-tree", "annotate"],
                             ["phylogenetic-tree", "delete"],
                             ["process", "read"],
                             ["process", "execute"]])


def initialize_database(flask_app):
    recreate_db = False
    if 'DB_CONNECTION_STRING' in flask_app.config:
        db_connection_string = flask_app.config['DB_CONNECTION_STRING']
        print(f"Connecting to {app_acronym.upper()} database server")
        print(db_connection_string)
        print("-----------------------------")
        if db_connection_string.startswith("sqlite://"):
            biobarcoding.engine = sqlalchemy.create_engine(db_connection_string,
                                                           echo=True,
                                                           connect_args={'check_same_thread': False},
                                                           poolclass=StaticPool)
        else:
            biobarcoding.engine = create_pg_database_engine(db_connection_string, app_acronym, recreate_db=recreate_db)

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


def initialize_postgis(flask_app):
    recreate_db = False
    if 'POSTGIS_CONNECTION_STRING' in flask_app.config:
        db_connection_string = flask_app.config['POSTGIS_CONNECTION_STRING']
        print(f"Connecting to {app_acronym}_geoserver database server")
        print(db_connection_string)
        print("-----------------------------")

        biobarcoding.postgis_engine = create_pg_database_engine(db_connection_string, f"{app_acronym}_geoserver",
                                                                recreate_db=recreate_db)
        # global DBSession # global DBSession registry to get the scoped_session
        DBSessionGeo.configure(
            bind=biobarcoding.postgis_engine)  # reconfigure the sessionmaker used by this scoped_session
        orm.configure_mappers()  # Important for SQLAlchemy-Continuum
        connection = biobarcoding.postgis_engine.connect()
        try:
            connection.execute("CREATE EXTENSION postgis")
        except:
            pass
        connection.execute("commit")
        connection.close()
        tables = ORMBaseGeo.metadata.tables
        connection = biobarcoding.postgis_engine.connect()
        table_existence = [biobarcoding.postgis_engine.dialect.has_table(connection, tables[t].name) for t in tables]
        connection.close()
        if False in table_existence:
            ORMBaseGeo.metadata.bind = biobarcoding.postgis_engine
            ORMBaseGeo.metadata.create_all()
        # Load base tables
        # initialize_gnd_geoserver_data() # not implemented
    else:
        print("No database connection defined (POSTGIS_CONNECTION_STRING), exiting now!")
        sys.exit(1)


def register_api(bp: Blueprint, view, endpoint: str, url: str, pk='id', pk_type='int'):
    view_func = view.as_view(endpoint)
    bp.add_url_rule(url, defaults={pk: None}, view_func=view_func, methods=['GET'])
    bp.add_url_rule(url, view_func=view_func, methods=['POST'])
    bp.add_url_rule(f'{url}<{pk_type}:{pk}>', view_func=view_func, methods=['GET', 'PUT', 'DELETE'])
    return view_func


def make_simple_rest_crud(entity, entity_name: str, execution_rules: Dict[str, str] = {},
                          aux_filter=None, default_filter=None, alt_getters: Dict[str, Callable] = None):
    """
    Create a CRUD RESTful endpoint

    :param entity: the entity class to manage
    :param entity_name: the name of the entity, in plural
    :param execution_rules: a dictionary of execution rules (according to the decorator "n_session") by CRUD method
    :param aux_filter: a function to augment the filtering capabilities (applies to GET without "id_")
    :param default_filter: A string with a default "filter" condition (applies to GET without "id_")
    :param alt_getters: a dictionary of alternative internal getters
    :return: the blueprint (to register it) and the class (if needed)
    """

    class CrudAPI(MethodView):

        @n_session(read_only=True, authr=execution_rules.get("r"))
        def get(self, _id=None):  # List or Read
            db = g.n_session.db_session
            r = ResponseObject()
            if _id is None:
                # List of all
                kwargs = parse_request_params()
                from ..services import get_query
                query, count = get_query(db, entity, aux_filter=aux_filter, default_filter=default_filter, **kwargs)
                # TODO Detail of fields
                r.count = count
                r.content = query.all()
            else:
                # Detail
                # TODO Detail of fields
                r.content = db.query(entity).filter(entity.id == int(_id)).first()
                r.count = 1

            return r.get_response()

        @n_session(authr=execution_rules.get("c"))
        def post(self):  # Create
            db = g.n_session.db_session
            r = ResponseObject()
            t = request.json
            s = entity.Schema().load(t, instance=entity(), partial=True)
            db.add(s)
            db.flush()
            r.content = s
            return r.get_response()

        @n_session(authr=execution_rules.get("u"))
        def put(self, _id):  # Update (total or partial)
            db = g.n_session.db_session
            r = ResponseObject()
            t = request.json
            if alt_getters and "put" in alt_getters:
                s = alt_getters["put"](db, _id)
            else:
                s = db.query(entity).filter(entity.id == int(_id)).first()
            s = entity.Schema().load(t, instance=s, partial=True)
            db.add(s)
            db.flush()
            r.content = s
            return r.get_response()

        @n_session(authr=execution_rules.get("d"))
        def delete(self, _id):  # Delete
            db = g.n_session.db_session
            r = ResponseObject()
            # TODO Multiple delete
            # if _id is None:
            #     # List of all
            #     query = db.query(entity)
            #     # Filter, Order, Pagination
            #     kwargs = check_request_params()
            #     if 'filter' in kwargs:
            #         query = query.filter(filter_parse(entity, kwargs.get('filter'))).delete()
            # else:
            if alt_getters and "delete" in alt_getters:
                s = alt_getters["delete"](db, _id)
            else:
                s = db.query(entity).filter(entity.id == int(_id)).first()
            db.delete(s)
            return r.get_response()

    # If the following line is uncommented, it complains on "overwriting an existing endpoint function".
    # This function is public, because it is just the schema, so, no problem.
    # @n_session()
    def get_entities_json_schema():
        r = ResponseObject()
        factory = SchemaFactory(StructuralWalker)
        r.content = factory(entity)

        return r.get_response()

    @n_session()
    def get_entity_json_schema(_id):
        db = g.n_session.db_session
        r = ResponseObject()
        factory = SchemaFactory(StructuralWalker)
        ent = db.query(entity).filter(entity.id == _id).first()
        r.content = dict(schema=factory(entity), data=ent)
        return r.get_response()

    # make_simple_rest_crud(entity, entity_name: str, execution_rules: Dict[str, str] = {})
    bp_entity = Blueprint(f'bp_{entity_name}', __name__)
    register_api(bp_entity, CrudAPI, entity_name, f"{app_api_base}/{entity_name}/", "_id", "string")
    # Two equivalent endpoints to obtain the JSON schema of an entity !!
    # (useful to automatically generate a Formly form)
    bp_entity.add_url_rule(f"{app_api_base}/{entity_name}.schema.json", view_func=get_entities_json_schema,
                           methods=['GET'])
    bp_entity.add_url_rule(f"{app_api_base}/{entity_name}/<int:_id>.schema.json", view_func=get_entity_json_schema,
                           methods=['GET'])

    return bp_entity, CrudAPI


# GENERIC REST FUNCTIONS

def chew_data(param):
    try:
        param = unquote(param)
    except Exception as e:
        pass
    try:
        param = json.loads(param)
    except Exception as e:
        pass
    return param


def decode_request_params(data):
    res = {}
    params = chew_data(data)
    for key in params:
        res[key] = chew_data(params[key])
    return res


def parse_request_params(data=None):
    kwargs = {'filter': [], 'order': [], 'pagination': {}, 'values': {}, 'searchValue': ''}
    if not data and request:
        if request.json:
            kwargs.update(parse_request_params(request.json))
        if request.data:
            kwargs.update(parse_request_params(request.data))
        if request.values:
            kwargs.update(parse_request_params(request.values))
    elif data:
        print(f'RAW_DATA: {data}')
        input = decode_request_params(data)
        for key in ('filter', 'order', 'pagination', 'values', 'searchValue'):
            try:
                i = input.pop(key)
            except Exception as e:
                continue
            kwargs[key] = i if i else kwargs[key]
        kwargs['values'].update(input)
    print(f'CLEAN_DATA: {kwargs}')
    return kwargs


def related_authr_ids(identity_id: int):
    # Get the organizations, groups, and roles associated to an identity
    ids = DBSession.query(Authorizable.id) \
        .join(OrganizationIdentity, OrganizationIdentity.organization_id==Authorizable.id, isouter=True) \
        .join(GroupIdentity, GroupIdentity.group_id==Authorizable.id, isouter=True) \
        .join(RoleIdentity, RoleIdentity.role_id==Authorizable.id, isouter=True) \
        .filter(or_(Authorizable.id==identity_id,
                    OrganizationIdentity.identity_id==identity_id,
                    GroupIdentity.identity_id==identity_id,
                    RoleIdentity.identity_id==identity_id))
    return [i for i, in ids]


def related_perm_ids(permission_id: int):
    # Get the mayor permissions that also allow permission_id
    ids = DBSession.query(PermissionType.id) \
        .filter(or_(PermissionType.id==permission_id,
                    PermissionType.rank <   # TODO: should be <= ?
                    DBSession.query(PermissionType.rank).filter(PermissionType.id==permission_id)))
    return [i for i, in ids]


def auth_filter(orm, permission_types_ids, object_types_ids,
                identity_id=None,
                object_uuids=None, time=None,
                permission_flag=False, authorizable_flag=False):
    """
    * orm: base (with uuid) to build the filter
    * identity_ids: who is requesting
    * permission_types_ids: what is requesting
    * object_types_ids: where is requesting
    * object_uuid: about what is requesting
    * time: when is requesting
    * permission_flag: extra info required
    * authorizable_flag: extra info required

    CollectionDetail (cd) <> Collection (c) > ACL <> ACLDetail (ad)

        SELECT
          CASE cd.object_uuid WHEN NULL THEN acl.object_uuid ELSE cd.object_uuid as oid,
          ad.permission_id as pid,
          ad.authr_id as aid  # To explain why it was authorized
        FROM CollectionDetail cd JOIN Collection c ON cd.collection_id=c.id RIGHT JOIN ACL ON c.uuid=acl.object_uuid JOIN ACLDetail ad ON acl.id=ad.acl_id
        WHERE ad.authorizable IN (<identities>)  # Authorizable
              AND object_type IN (<collection object type>[, <target object type>])  # Object types
              AND date time BETWEEN ad.validity_start AND ad.validity_end
              [AND permission_types IN (...)]
              [AND object_uuids IN (...)]

        * Añadir lista object_ids objeto sirve para ver permisos de esos objetos concretos
        * Añadir lista id tipos permisos sirve para ver si esos tipos están

    @return: <orm_clause_filter> || object_uuids[, permissions][, authorizables]
    """
    from flask import current_app
    if current_app.config["ACL_ENABLED"] != 'True':
        return True

    from datetime import datetime
    time = time if time else datetime.now()
    filter_clause = [
        or_(time >= ACLDetail.validity_start, ACLDetail.validity_start == None),
        or_(time <= ACLDetail.validity_end, ACLDetail.validity_end == None)
    ]
    filter_clause.append(ACL.object_type.in_(object_types_ids)) if object_types_ids else None
    filter_clause.append(ACL.object_uuid.in_(object_uuids)) if object_uuids else None
    try:
        # By identity or authorizables associated with the identity
        if isinstance(identity_id, int):
            filter_clause.append(ACLDetail.authorizable_id.in_(related_authr_ids(identity_id)))
        elif isinstance(identity_id, (tuple, list, set)):
            filter_clause.append(ACLDetail.authorizable_id.in_(identity_id))
        else:
            filter_clause.append(ACLDetail.authorizable_id.in_(related_authr_ids(g.n_session.identity.id)))
        # By permission or superior permissions
        if isinstance(permission_types_ids, int):
            filter_clause.append(ACLDetail.permission_id.in_(related_perm_ids(permission_types_ids)))
        elif isinstance(permission_types_ids, (tuple, list, set)):
            filter_clause.append(ACLDetail.permission_id.in_(permission_types_ids))

        collected = DBSession.query(CollectionDetail.object_uuid) \
            .join(Collection) \
            .join(ACL, Collection.uuid==ACL.object_uuid) \
            .join(ACLDetail) \
            .filter(*filter_clause)
        uncollected = DBSession.query(ACL.object_uuid).join(ACLDetail).filter(*filter_clause)

        final_query = collected.union(uncollected)
        if permission_flag or authorizable_flag:
            entities = []
            entities += [ACLDetail.permission] if permission_flag else []
            entities += [ACLDetail.authorizable] if authorizable_flag else []
            return final_query.with_entities(entities)
        else:
            return orm.uuid.in_(final_query)
    except Exception as e:
        print(e)
        return None


def filter_parse(orm, filter, aux_filter=None, session=None):
    """
    @param orm: <ORMSqlAlchemy>
    @param filter:
        [{"field1": "<condition>", "field2": "<condition>"}, {"field1": "<condition"}]
        <condition>: {"op": "<operador", "left": "<valor>", "right": <valor>, "unary": <valor>}
    @param aux_filter: <callable function(filter)>
    @return: <orm_clause_filter>
    """

    def get_condition(orm, field, condition):
        obj = getattr(orm, field)
        op = condition["op"]
        value, left, right = condition.get("unary"), condition.get("left"), condition.get("right")
        if op == "out":
            return obj.notin_(value)
        if op == "in":
            return obj.in_(value)
        if op == "eq":
            return obj == value
        if op == "between":
            return obj.between_(left, right)
        if op == "le":
            return obj <= value
        if op == "ge":
            return obj >= value
        if op == "contains":
            return value.in_(obj)

        return True

    try:
        if not isinstance(filter, (list, tuple)):
            filter = [filter]
        or_clause = []
        for clause in filter:
            and_clause = []
            for field, condition in clause.items():
                if hasattr(orm, field):
                    and_clause.append(get_condition(orm, field, condition))
                else:
                    print(f'Warning: Unknown column "{field}" for the given tables.')
            if aux_filter:
                try:
                    lst = aux_filter(clause, session)
                except:
                    lst = aux_filter(clause)
                and_clause += lst

            or_clause.append(and_(*and_clause))
        return or_(*or_clause)
    except Exception as e:
        traceback.print_exc()
        return None


def order_parse(orm, sort, aux_order=None):
    """
    @param orm: <ORMSqlAlchemy>
    @param sort: [{"order": "<asc|desc>", "field": "<name>"}]
    @param aux_order: <callable function(order)>
    @return: <orm_clause_order>
    """

    def get_condition(orm, clause):
        obj = getattr(orm, clause["field"])
        if clause["order"].startswith("asc"):
            return obj.asc()
        if clause["order"].startswith("desc"):
            return obj.desc()
        return True

    try:
        if not isinstance(sort, (list, tuple)):
            order = [sort]
        else:
            order = sort
        ord_clause = []
        for clause in order:
            if hasattr(orm, clause.get("field", "")):
                ord_clause.append(get_condition(orm, clause))
            else:
                print(f'Unknown column "{clause.get("field", "")}" for the given tables.')
                if aux_order:
                    ord_clause += aux_order(clause)
        return ord_clause
    except Exception as e:
        print(e)
        return None


def get_content(session, clazz, issues, id_=None, aux_filter=None, query=None):
    content = None
    count = 0
    try:
        if id_ is None:
            # List of all
            kwargs = parse_request_params()
            from ..services import get_query
            query2, count = get_query(session, clazz, query, aux_filter=aux_filter, **kwargs)
            content = query2.all()
        else:
            # Detail
            if not query:
                query = session.query(clazz)
            content = query.filter(clazz.id == id_).first()
    except Exception as e:
        traceback.print_exc(e)
        _, status = issues.append(Issue(IType.ERROR,
                                        f'READ "{clazz.__name__}" data: The "{clazz.__name__}"'
                                        f' data could not be read.')), 500
    if not content:
        _, status = issues.append(Issue(IType.INFO, f'no data available')), 200
    else:
        _, status = issues.append(
            Issue(IType.INFO, f'READ "{clazz.__name__}": The "{clazz.__name__}"'
                              f' data were successfully read')), 200
    return issues, content, count, status


def initialize_ssh(flask_app):
    """
    ssh_resources = conn.execute("SELECT * FROM jobs_compute_resources join jobs_job_mgmt_types as jt" +
                                 " on jm_type_id=jt.id where jt.name='ssh'").fetchall()
    :param flask_app:
    :return:
    """
    session = DBSession()
    ssh_resources = session.query(ComputeResource, JobManagementType).filter(
        ComputeResource.jm_type_id == JobManagementType.id, JobManagementType.name == 'ssh').all()
    for ssh_res in ssh_resources:
        compute_resource = ssh_res.ComputeResource
        ssh_credentials = compute_resource.jm_credentials
        ssh_location = compute_resource.jm_location
        cmd = ["ssh", "-o", "StrictHostKeyChecking=no",
               "-o", f"UserKnownHostsFile={ssh_credentials['known_hosts_filepath']}",
               f"{ssh_credentials['username']}@{ssh_location['host']}", "'echo'"]
        subprocess.run(cmd)


# GALAXY INIZIALIZATION
def install_galaxy_tools(gi, tools):
    '''
    tool_shed_url, name, owner,
    changeset_revision,
    install_tool_dependencies = False,
    install_repository_dependencies = False,
    install_resolver_dependencies = False,
    tool_panel_section_id = None,
    new_tool_panel_section_label = Non
    '''
    if isinstance(tools, list):
        for tool in tools:
            tool_shed_url = 'https://' + tool['tool_shed']
            gi.toolshed.install_repository_revision(tool_shed_url=tool_shed_url,
                                                    name=tool['name'],
                                                    owner=tool['owner'],
                                                    changeset_revision=tool['changeset_revision'],
                                                    install_tool_dependencies=True,
                                                    install_repository_dependencies=True,
                                                    install_resolver_dependencies=True,
                                                    new_tool_panel_section_label='New'
                                                    )


def check_galaxy_tools(wf1_dic, wf2_dic):
    steps1 = wf1_dic['steps']
    steps2 = wf2_dic['steps']
    tool_list = list()
    for step, content in steps1.items():
        if 'errors' in content:
            if content['errors'] == "Tool is not installed":
                # TODO depende de la versión de galaxy esto lleva un punto al final o no xq lo que hay que buscar
                #  otra cosa
                tool_list.append(steps2[step]['tool_shed_repository'])
    if len(tool_list) == 0:
        return 'all tools are installed'
    else:
        return tool_list


def initialize_galaxy(flask_app):
    session = DBSession()
    galaxy_resources = session.query(ComputeResource, JobManagementType).filter(
        ComputeResource.jm_type_id == JobManagementType.id, JobManagementType.name == 'galaxy').all()
    for galaxy_res in galaxy_resources:
        compute_resource = galaxy_res.ComputeResource
        url = compute_resource.jm_location['url']
        api_key = compute_resource.jm_credentials['api_key']

        # Install basic workflow if it is not installed
        gi = galaxy.GalaxyInstance(url=url, key=api_key)
        path = ROOT + '/biobarcoding/workflows/'
        for workflow in os.listdir(path):
            workflow_path = path + workflow
            with open(workflow_path, 'r') as f:
                wf_dict_in = json.load(f)
            name = wf_dict_in['name']
            wf = gi.workflows.get_workflows(name=name)
            w_id = wf[0].get('id') if wf else None
            if w_id:
                gi.workflows.delete_workflow(w_id)
            wf = gi.workflows.import_workflow_from_local_path(workflow_path)
            wf_dict_out = gi.workflows.export_workflow_dict(wf['id'])
            list_of_tools = check_galaxy_tools(wf_dict_out, wf_dict_in)
            install_galaxy_tools(gi, list_of_tools)
        # conversion of workflows galaxy into workflows formly
        # TODO convert_workflows_to_formly()


def initialize_database_chado(flask_app):
    cfg = flask_app.config
    if 'CHADO_CONNECTION_STRING' in flask_app.config:
        db_connection_string = cfg["CHADO_CONNECTION_STRING"]
    elif 'CHADO_DATABASE' in flask_app.config:
        db_connection_string = f'postgresql://{cfg["CHADO_USER"]}:{cfg["CHADO_PASSWORD"]}@{cfg["CHADO_HOST"]}:{cfg["CHADO_PORT"]}/{cfg["CHADO_DATABASE"]}'
        print("Connecting to Chado database server")
        print(db_connection_string)
        print("-----------------------------")
        biobarcoding.chado_engine = sqlalchemy.create_engine(db_connection_string, echo=False)
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


def initialize_chado_edam(flask_app):
    from sqlalchemy import text, Table, Index, MetaData
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

        db_connection_string = f'postgresql://{cfg["CHADO_USER"]}:{cfg["CHADO_PASSWORD"]}@{cfg["CHADO_HOST"]}:{cfg["CHADO_PORT"]}/{cfg["CHADO_DATABASE"]}'
        print("Connecting to Chado database server to insert EDAM")
        print(db_connection_string)
        print("-----------------------------")
        biobarcoding.chado_engine = sqlalchemy.create_engine(db_connection_string, echo=True)

        with biobarcoding.chado_engine.connect() as conn:
            relationship_id = conn.execute(
                text("select * from INFORMATION_SCHEMA.TABLES where TABLE_NAME = 'db_relationship'")).fetchone()
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
                                                   ForeignKey('cvterm.cvterm_id', ondelete="CASCADE",
                                                              initially="DEFERRED"),
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
                    conn.execute(
                        text("delete from cvterm where name = 'has_function'and cv_id = 4 and is_obsolete = 0"))
                    conn.execute(text("delete from cvterm where name = 'has_input'and cv_id = 4 and is_obsolete = 0"))
                    conn.execute(text("delete from cvterm where name = 'has_output'and cv_id = 4 and is_obsolete = 0"))
                    conn.execute(text("delete from cvterm where name = 'is_a'and cv_id = 4 and is_obsolete = 0"))
                    conn.execute(
                        text("delete from cvterm where name = 'is_anonymous'and cv_id = 6 and is_obsolete = 0"))

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
                    # We commit 2 times because we need a commit on the last block to execute this one
                    with biobarcoding.chado_engine.connect() as conn:
                        # FILL THE db_relationship TABLE
                        edam_id = conn.execute(text("select db_id from db where name=\'EDAM\'")).fetchone()[0]
                        db_ids = conn.execute(text("select db_id from db where db_id > " + str(edam_id))).fetchall()
                        type_id = conn.execute(text(
                            "select cvterm_id from cvterm join cv on cvterm.cv_id = cv.cv_id where cvterm.name = \'part_of\' "
                            "and cv.name = \'relationship\'")).fetchone()[0]

                        for db_id in db_ids:
                            conn.execute(text("INSERT INTO db_relationship (type_id,subject_id,object_id)" +
                                              "VALUES(" + str(type_id) + "," + str(db_id[0]) + "," + str(
                                edam_id) + ")"))
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
