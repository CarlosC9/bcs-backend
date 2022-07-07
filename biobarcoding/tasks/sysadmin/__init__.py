from ..system import SA_TASK_SESSION
from ... import get_global_configuration_variable
from ...rest import app_api_base
from ...services import log_exception

ENDPOINT = get_global_configuration_variable("ENDPOINT_URL")
COOKIES_FILE_PATH = get_global_configuration_variable("COOKIES_FILE_PATH")
REQUEST_URL = f"{ENDPOINT}{app_api_base}"


def create_request(url_suffix, **kwargs):
    try:
        url = f"{REQUEST_URL}{url_suffix}"
        print('POST ' + url)
        return SA_TASK_SESSION.post(url, json=kwargs, headers={'Content-Type': 'application/json'})
    except Exception as e:
        print(f'Something went wrong when creating the {kwargs.get("name")}. It may already exist.')
        log_exception(e)
        return None


def read_request(url_suffix, **kwargs):
    try:
        url = f"{REQUEST_URL}{url_suffix}"
        print('GET ' + url)
        return SA_TASK_SESSION.get(url, json=kwargs)
    except Exception as e:
        print(f'Something went wrong when reading the {kwargs.get("name")}.')
        log_exception(e)
        return None
