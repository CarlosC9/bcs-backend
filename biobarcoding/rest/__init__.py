import json
import logging
import os
import sys
from enum import Enum
from typing import Dict, List
from urllib.parse import unquote

import redis
from alchemyjsonschema import SchemaFactory, StructuralWalker
from attr import attrs, attrib
from bioblend import galaxy
from flask import Response, Blueprint, g, request
from flask.views import MethodView
from sqlalchemy import orm, and_, or_
from sqlalchemy.pool import StaticPool

import biobarcoding
from ..authentication import n_session
from ..common import generate_json, ROOT
from ..common.helpers import get_module_logger
from ..common.pg_helpers import create_pg_database_engine, load_table, load_many_to_many_table, \
    load_computing_resources, load_processes_in_computing_resources, load_process_input_schema, load_table_extended
from ..db_models import DBSession, DBSessionChado, ORMBaseChado, DBSessionGeo
from ..db_models.bioinformatics import *
from ..db_models.geographics import *
from ..db_models.jobs import *
from ..db_models.sysadmin import *
from ..rest.socket_service import SocketService

app_api_base = "/api"  # Base for all RESTful calls
app_gui_base = "/gui"  # Base for the Angular2 GUI
app_external_gui_base = "/gui_external"  # Base for the Angular2 GUI when loaded from URL
app_proxy_base = "/pxy"  # Base for the BCS Proxy

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
    app_dirs = AppDirs("bcs-backend")
    # Default directories, multi-platform
    data_path = app_dirs.user_data_dir
    cache_path = app_dirs.user_cache_dir

    REDIS_HOST = "redis"
    REDIS_PORT = 6379
    BROKER_URL = f"redis://{REDIS_HOST}:{REDIS_PORT}/0"
    BACKEND_URL = BROKER_URL

    return dict(
        # BCS (SYSTEM DB)
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
        # GEOSERVER address from BCS-Backend
        GEOSERVER_URL="localhost:9180",
        # GEO (GEOSPATIAL DATA)
        GEOSERVER_USER="admin",
        GEOSERVER_PASSWORD="ngd_ad37",
        GEOSERVER_HOST="geoserver",
        GEOSERVER_PORT="8080",
        # PostGIS address from BCS-Backend
        POSTGIS_CONNECTION_STRING="postgresql://postgres:postgres@localhost:5435/",
        # PostGIS address from Geoserver ("host", "port" could differ from those in POSTGIS_CONNECTION_STRING)
        POSTGIS_USER="postgres",
        POSTGIS_PASSWORD="postgres",
        POSTGIS_PORT="5432",
        POSTGIS_HOST="localhost",
        POSTGIS_DB="ngd_geoserver",
        # COMPUTE RESOURCES
        RESOURCES_CONFIG_FILE_PATH="/home/resources_config.json",
        JOBS_LOCAL_WORKSPACE=os.path.expanduser('~/ngd_jobs'),
        SSH_JOBS_DEFAULT_REMOTE_WORKSPACE="/tmp",
        GALAXY_API_KEY="fakekey",
        GALAXY_LOCATION="http://localhost:8080",
        # MISC
        GOOGLE_APPLICATION_CREDENTIALS=f"{data_path}/firebase-key.json",  # Firebase
        ENDPOINT_URL="http://localhost:5000",  # Self "bcs-backend" address (so Celery can update things)
        COOKIES_FILE_PATH="/tmp/bcs-cookies.txt",  # Where cookies are stored by Celery
        CACHE_FILE_LOCATION=f"{cache_path}/cache",  # Cached things
        SAMESITE_NONE="True",
        TESTING="True",
        SELF_SCHEMA="",
    )


def prepare_default_configuration(create_directories):
    default_cfg = get_default_configuration_dict()

    from appdirs import AppDirs
    app_dirs = AppDirs("bcs-backend")

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
           data_path + os.sep + "bcs_local.conf"


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
                    print(f"WARNING: {key.strip()} is not in the default config")
        for key in default_cfg_keys:
            print(f"{key} inserted to {file_name} with value {default_cfg[key]}")
            file.write(f'\n{key}="{default_cfg[key]}"')


def load_configuration_file(flask_app):
    # Initialize configuration
    try:
        _, file_name = prepare_default_configuration(False)
        if not os.environ.get(biobarcoding.config_file_var):
            logger.debug(f"Trying {file_name}")
            if os.path.isfile(file_name):
                print(f"Assuming {file_name} as configuration file")
                logger.debug(f"Assuming {file_name} as configuration file")
                found = True
                complete_configuration_file(file_name)
                os.environ[biobarcoding.config_file_var] = file_name
            else:
                found = False
                logger.debug(f"Creating {file_name} as configuration file")
            if not found:
                cfg, file_name = prepare_default_configuration(True)
                print(f"Generating {file_name} as configuration file:\n{cfg}")
                with open(file_name, "wt") as f:
                    f.write(cfg)
                os.environ[biobarcoding.config_file_var] = file_name
        else:
            complete_configuration_file(os.environ.get(biobarcoding.config_file_var))
        print("-----------------------------------------------")
        print(f'Configuration file at: {os.environ[biobarcoding.config_file_var]}')
        print("-----------------------------------------------")
        flask_app.config.from_envvar(biobarcoding.config_file_var)
        logger.debug(flask_app.config)

    except Exception as e:
        print(f"{biobarcoding.config_file_var} environment variable not defined!")
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


# SOCKET INITIALIZATION
def init_socket(socketio):
    SocketService(socketio)


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

tm_processes = {  # Preloaded processes
    "c8df0c20-9cd5-499b-92d4-5fb35b5a369a": "MSA ClustalW",
    "ec40143f-ae32-4dac-9cfb-caa047e1adb1": "ClustalW-PhyMl",
    "c87f58b6-cb06-4d39-a0b3-72c2705c5ae1": "PAUP Parsimony",
    "3e0240e8-b978-48a2-8fdd-9f31f4264064": "Phylogenetic Diversity Analyzer",
    "c55280d0-f916-4401-a1a4-bb26d8179fd7": "MSA ClustalW + PAUP Parsimony",
    "ce018826-7b20-4b70-b9b3-168c0ba46eec": "PAUP Parsimony + Phylogenetic Diversity Analyzer",
    "5b315dc5-ad12-4214-bb6a-bf013f0e4b8c": "MSA ClustalW + PAUP Parsimony + Phylogenetic Diversity Analyzer",
}

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

tm_browser_filter_form_fields = ("id", "uuid")
tm_browser_filter_forms = [
    (bio_object_type_id["sequence"], "6bb4cbbe-d6cd-4c2e-bb59-782e3c9e9f6c"),
    (bio_object_type_id["multiple-sequence-alignment"], "44363784-d304-4e7e-a507-8bae48598e50"),
    (bio_object_type_id["phylogenetic-tree"], "8b62f4aa-d32a-4841-89f5-9ed50da44121"),
]

# 02f44e54-f139-4ea0-a1bf-fe27054c0d6c
# 903a73a9-5a4e-4cec-b8fa-4fc9bd5ffab5
# 5c4ba6db-e7f2-4d5c-a89a-76059ac116b1
# ea647c4e-2063-4246-bd9a-42f6a57fb9ea
# 985c01ca-d9d2-4df5-a8b9-8a6da251d7d4
# f167eac0-2a23-4e74-bb1c-abdfb5f74a92
# 4cfcd389-ed9e-4174-aa99-150f176e8eec
# caaca280-2290-4625-b5c0-76bcfb06e9ac
# 15aa399f-dd58-433f-8e94-5b2222cd06c9
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
    # load_table_extended(DBSession, BrowserFilterForm, tm_browser_filter_form_fields, tm_browser_filter_forms)

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
    load_many_to_many_table(DBSession, ObjectTypePermissionType, ObjectType, PermissionType,
                            ["object_type_id", "permission_type_id"],
                            [["sequence", "read"], ["sequence", "annotate"], ["sequence", "delete"],
                             ["multiple-sequence-alignment", "read"], ["multiple-sequence-alignment", "annotate"],
                             ["multiple-sequence-alignment", "delete"],
                             ["phylogenetic-tree", "read"], ["phylogenetic-tree", "annotate"],
                             ["phylogenetic-tree", "delete"],
                             ["process", "read"], ["process", "execute"]])

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


def initialize_postgis(flask_app):
    recreate_db = False
    if 'POSTGIS_CONNECTION_STRING' in flask_app.config:
        db_connection_string = flask_app.config['POSTGIS_CONNECTION_STRING']
        print("Connecting to ngd_geoserver database server")
        print(db_connection_string)
        print("-----------------------------")

        biobarcoding.postgis_engine = create_pg_database_engine(db_connection_string, "ngd_geoserver",
                                                                recreate_db=recreate_db)
        # global DBSession # global DBSession registry to get the scoped_session
        DBSessionGeo.configure(
            bind=biobarcoding.postgis_engine)  # reconfigure the sessionmaker used by this scoped_session
        orm.configure_mappers()  # Important for SQLAlchemy-Continuum
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


def make_simple_rest_crud(entity, entity_name: str, execution_rules: Dict[str, str] = {}):
    """
    Create a CRUD RESTful endpoint

    :param entity: the entity class to manage
    :param entity_name: the name of the entity, in plural
    :param execution_rules: a dictionary of execution rules (according to the decorator "bcs_session") by CRUD method
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
                query, count = get_query(db, entity, **kwargs)
                # TODO Detail of fields
                r.content = query.all()
            else:
                # Detail
                # TODO Detail of fields
                r.content = db.query(entity).filter(entity.id == _id).first()

            return r.get_response()

        @n_session(authr=execution_rules.get("c"))
        def post(self):  # Create
            db = g.n_session.db_session
            r = ResponseObject()
            t = request.json
            s = entity.Schema().load(t, instance=entity(), partial=True)
            db.add(s)
            return r.get_response()

        @n_session(authr=execution_rules.get("u"))
        def put(self, _id):  # Update (total or partial)
            db = g.n_session.db_session
            r = ResponseObject()
            t = request.json
            s = db.query(entity).filter(entity.id == _id).first()
            s = entity.Schema().load(t, instance=s, partial=True)
            db.add(s)
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
            s = db.query(entity).filter(entity.id == _id).first()
            db.delete(s)
            return r.get_response()

    # If the following line is uncommented, it complains on "overwriting an existing endpoint function". This function is public, because it is just the schema, so, no problem.
    # @bcs_session()
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

    bp_entity = Blueprint(f'bp_{entity_name}', __name__)
    register_api(bp_entity, CrudAPI, entity_name, f"{app_api_base}/{entity_name}/", "_id")
    bp_entity.add_url_rule(f"{app_api_base}/{entity_name}.schema.json", view_func=get_entities_json_schema,
                           methods=['GET'])
    bp_entity.add_url_rule(f"{app_api_base}/{entity_name}/<int:_id>.schema.json", view_func=get_entity_json_schema,
                           methods=['GET'])

    return bp_entity, CrudAPI


# GENERIC REST FUNCTIONS

def decode_request_params(data):
    res = {}
    for key in data:
        value = data[key]
        try:
            value = unquote(value)
        except Exception as e:
            pass
        try:
            value = json.loads(value)
        except Exception as e:
            pass
        res[key] = value
    return res


def parse_request_params(data=None):
    kwargs = {'filter': [], 'order': [], 'pagination': {}, 'value': {}, 'searchValue': ''}
    if not data:
        if request.json:
            kwargs.update(parse_request_params(request.json))
        if request.values:
            kwargs.update(parse_request_params(request.values))
    else:
        print(f'DATA: {data}')
        input = decode_request_params(data)
        for key in ('filter', 'order', 'pagination', 'value', 'searchValue'):
            try:
                i = input.pop(key)
            except Exception as e:
                continue
            kwargs[key] = i if i else kwargs[key]
        kwargs['value'].update(input)
    print(f'KWARGS: {kwargs}')
    return kwargs


def related_authr_ids(engine, identity_id) -> List[int]:
    pass


def auth_filter():
    """
    * Filter
    * Identity (ids) 
    * object types (ids)
    * permission types (ids)
    * object ids
    * date time
    
    CollectionDetail (cd) <> Collection (c) > ACL <> ACLDetail (ad)
    
SELECT 
  CASE cd.object_id WHEN NULL THEN acl.object_uuid ELSE cd.object_uuid as oid, 
  ad.permission_id as pid, 
  ad.authr_id as aid  # To explain why it was authorized
FROM CollectionDetail cd JOIN Collection c ON cd.col_id=c.id RIGHT JOIN ACL ON c.uuid=acl.object_id JOIN ACLDetail ad ON acl.id=ad.acl_id
WHERE ad.authorizable IN (<identities>)  # Authorizable 
      AND object_type IN (<collection object type>[, <target object type>])  # Object types
      AND date time BETWEEN ad.validity_start AND ad.validity_end
      [AND permission_types IN (...)]
      [AND object_uuids IN (...)]
      
* Añadir lista ids objeto sirve para ver permisos de esos objetos concretos
* Añadir lista id tipos permisos sirve para ver si esos tipos están
    
    @return: 
    """
    pass


def filter_parse(orm, filter, aux_filter=None):
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
                    print(f'Unknown column "{field}" for the given tables.')
            if aux_filter:
                and_clause += aux_filter(clause)
            or_clause.append(and_(*and_clause))
        return or_(*or_clause)
    except Exception as e:
        print(e)
        return None


def order_parse(orm, sort, aux_order=None):
    """
    @param orm: <ORMSqlAlchemy>
    @param sort: [{"order": "<asc|desc>", "field": "<name>"}]
    @param aux_order: <callable function(order)>
    @return: <orm_clause_order>
    """

    def get_condition(orm, clause):
        obj = getattr(orm, clause.field)
        if clause.order == "asc":
            return obj.asc()
        if clause.order == "desc":
            return obj.desc()
        return True

    try:
        if not isinstance(sort, (list, tuple)):
            order = [sort]
        ord_clause = []
        for clause in order:
            if hasattr(orm, clause.field):
                ord_clause.append(get_condition(orm, clause))
            else:
                print(f'Unknown column "{clause.field}" for the given tables.')
                if aux_order:
                    ord_clause += aux_order(clause)
        return ord_clause
    except Exception as e:
        print(e)
        return None
