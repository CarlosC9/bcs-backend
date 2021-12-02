# GLOBAL VARIABLES
import configparser
import os
from typing import Dict

flask_app = None  # Flask app
engine = None  # SQLAlchemy database engine
postgis_engine = None
chado_engine = None  # SQLAlchemy database engine
celery = None
case_sensitive = True  # Define whether comparison of user defined identifiers is case sensitive or not
config_file_var = "BCS_CONFIG_FILE"
global_configuration = None


def remove_quotes(s: str) -> str:
    return s.lstrip('\'"').rstrip('\'"')


def expand_paths(key: str, value: str) -> str:
    if key.endswith('_dir') or key.endswith('_location') or key.endswith('_file'):
        return os.path.expanduser(value)
    else:
        return value


def get_global_configuration_variable(key: str, default=None) -> str:
    def read_configuration() -> Dict[str, str]:
        """
        If environment variable "MAGIC_NIS_SERVICE_CONFIG_FILE" is defined, and the contents is the name of an existing file,
        read it as a configuration file and return the result

        :return:
        """
        if os.environ.get(config_file_var):
            fname = os.environ[config_file_var]
            if os.path.exists(fname):
                with open(fname, 'r') as f:
                    config_string = '[asection]\n' + f.read()
                config = configparser.ConfigParser()
                config.read_string(config_string)
                return {k: expand_paths(k, remove_quotes(v)) for k, v in config.items("asection")}
        return {}

    global global_configuration
    if global_configuration is None:
        global_configuration = read_configuration()
    if default is None:
        from .rest import get_default_configuration_dict
        default = get_default_configuration_dict().get(key)
    return global_configuration.get(key.lower(), default)


app_acronym = "bcs"  # Application Short acronym (internal use)