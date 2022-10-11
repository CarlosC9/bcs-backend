import json
from os import path

from ..system import SA_TASK_SESSION
from ... import get_global_configuration_variable
from ...rest import app_api_base
from ...services import log_exception

ENDPOINT = get_global_configuration_variable("ENDPOINT_URL")
COOKIES_FILE_PATH = get_global_configuration_variable("COOKIES_FILE_PATH")
REQUEST_URL = f"{ENDPOINT}{app_api_base}"


def load_response(r):
    try:
        return json.loads(r.text)
    except Exception as e:
        return {}


def get_response_count(r):
    return load_response(r).get('count')


def get_response_content(r):
    return load_response(r).get('content')


def get_response_id(r):
    try:
        c = get_response_content(r)
        if isinstance(c, (tuple, list, set)):
            c = c[0] if len(c) == 1 else {}
        return c.get('id')
    except Exception as e:
        print('missing id for ', c)
        return None


def create_request(url_suffix, **kwargs):
    try:
        url = path.join(REQUEST_URL, url_suffix)
        print('POST ' + url)
        return SA_TASK_SESSION.post(url, json=kwargs, headers={'Content-Type': 'application/json'})
    except Exception as e:
        print(f'Something went wrong when creating the {kwargs.get("name")}. It may already exist.')
        log_exception(e)
        return None


def read_request(url_suffix, **kwargs):
    try:
        url = path.join(REQUEST_URL, url_suffix)
        print('GET ' + url)
        return SA_TASK_SESSION.get(url, json=kwargs)
    except Exception as e:
        print(f'Something went wrong when reading the {kwargs.get("name")}.')
        log_exception(e)
        return None


def update_request(url_suffix, **kwargs):
    try:
        url = path.join(REQUEST_URL, url_suffix)
        print('PUT ' + url)
        return SA_TASK_SESSION.put(url, json=kwargs)
    except Exception as e:
        print(f'Something went wrong when reading the {kwargs.get("name")}.')
        log_exception(e)
        return None
