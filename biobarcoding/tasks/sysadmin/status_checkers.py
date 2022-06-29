import subprocess

from . import *
from ...services import log_exception


def run():
	try:
		url = f"{REQUEST_URL}/sys/status_checkers"
		cmd = ["curl", "--cookie", COOKIES_FILE_PATH, "--cookie-jar", COOKIES_FILE_PATH, "-X", "PUT", url]
		return subprocess.run(cmd).stdout
	except Exception as e:
		log_exception(e)
		return None
