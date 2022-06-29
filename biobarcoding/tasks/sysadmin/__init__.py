from ... import get_global_configuration_variable
from ...rest import app_api_base

ENDPOINT = get_global_configuration_variable("ENDPOINT_URL")
COOKIES_FILE_PATH = get_global_configuration_variable("COOKIES_FILE_PATH")
REQUEST_URL = f"{ENDPOINT}{app_api_base}"
