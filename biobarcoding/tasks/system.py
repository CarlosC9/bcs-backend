import subprocess

from . import celery_app
from .. import get_global_configuration_variable
from ..rest import app_api_base
from ..services import log_exception

SA_PROC_PATH = 'biobarcoding.tasks.sysadmin.'
CELERY_USER = "celery_user"
ENDPOINT = get_global_configuration_variable("ENDPOINT_URL")
COOKIES_FILE_PATH = get_global_configuration_variable("COOKIES_FILE_PATH")


# ----------------------------------------------------------------------------------------------------------------------
# THE TASKS
# ----------------------------------------------------------------------------------------------------------------------

@celery_app.on_after_finalize.connect
def setup_periodic_tasks(sender, **kwargs):
    sender.add_periodic_task(30.0, sa_task.s('status_checkers'), name='check subsystem status')


def api_login():
    url = f"{ENDPOINT}{app_api_base}/authn?user={CELERY_USER}"
    cmd = ["curl", "--cookie", COOKIES_FILE_PATH, "--cookie-jar", COOKIES_FILE_PATH, "-X", "PUT", url]
    subprocess.run(cmd)


def api_logout():
    url = f"{ENDPOINT}{app_api_base}/authn?user={CELERY_USER}"
    cmd = ["curl", "--cookie", COOKIES_FILE_PATH, "--cookie-jar", COOKIES_FILE_PATH, "-X", "DELETE", url]
    subprocess.run(cmd)


@celery_app.task(name="sa_task", acks_late=True)
def sa_task(process: str):
    """
    Launch a specific sysadmin task
    :param process:
    :return:
    """

    try:
        api_login()
        print('sa_task: ' + process)
        from importlib import import_module
        process_module = import_module(SA_PROC_PATH + process)
        _ = process_module.run()
        api_logout()
        return _
    except Exception as e:
        log_exception(e)
        return None
