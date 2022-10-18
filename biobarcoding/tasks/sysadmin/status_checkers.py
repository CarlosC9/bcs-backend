import json

from . import REQUEST_URL
from ..system import SA_TASK_SESSION, sa_task
from ...services import log_exception


def run():
	try:
		_ = json.loads(SA_TASK_SESSION.put(f"{REQUEST_URL}/sys/status_checkers").text)
		_.pop('content')
		print(_)
		return 'DONE'
	except Exception as e:
		log_exception(e)
		return 'EXCEPTION'
	finally:
		sa_task.apply_async(args=('status_checkers',), countdown=30)

