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

    def define_config_file_var():
        from biobarcoding.rest import complete_configuration_file, prepare_default_configuration
        _, file_name = prepare_default_configuration(False)
        if not os.environ.get(config_file_var):
            if os.path.isfile(file_name):
                found = True
                complete_configuration_file(file_name)
            else:
                found = False
            if not found:
                cfg, file_name = prepare_default_configuration(True)
                with open(file_name, "wt") as f:
                    f.write(cfg)
        else:
            file_name = os.environ.get(config_file_var)
            complete_configuration_file(file_name)
        return file_name

    def read_configuration() -> Dict[str, str]:
        """
        If environment variable "BCS_CONFIG_FILE" is defined, and the contents is the name of an existing file,
        read it as a configuration file and return the result

        :return:
        """
        if not os.environ.get(config_file_var):
            fname = define_config_file_var()
            os.environ[config_file_var] = fname
        else:
            fname = os.environ[config_file_var]
        if fname.strip() != "":
            if os.path.exists(fname):
                with open(fname, 'r') as f:
                    config_string = '[asection]\n' + f.read()
                config = configparser.ConfigParser()
                config.read_string(config_string)
                return {k: expand_paths(k, remove_quotes(v)) for k, v in config.items("asection")}
            else:
                return {}
        else:
            return {}

    global global_configuration
    if global_configuration is None:
        global_configuration = read_configuration()
    if default is None:
        from .rest import get_default_configuration_dict
        default = get_default_configuration_dict().get(key)
    return global_configuration.get(key.lower(), default)


app_acronym = "ngd"  # Application Short acronym (internal use)