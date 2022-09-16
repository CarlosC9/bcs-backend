import requests

from . import celery_app
from .. import get_global_configuration_variable
from ..rest import app_api_base
from ..services import log_exception

SA_PROC_PATH = 'biobarcoding.tasks.sysadmin.'
CELERY_USER = "celery_user"
ENDPOINT = get_global_configuration_variable("ENDPOINT_URL")

if 'SA_TASK_SESSION' not in globals():
    SA_TASK_SESSION = requests.Session()


# ----------------------------------------------------------------------------------------------------------------------
# THE TASKS
# ----------------------------------------------------------------------------------------------------------------------

@celery_app.on_after_finalize.connect
def setup_periodic_tasks(sender, **kwargs):
    sender.add_periodic_task(30.0, sa_task.s('status_checkers'), name='check subsystem status')


def sa_task_login():
    SA_TASK_SESSION.put(f"{ENDPOINT}{app_api_base}/authn?user={CELERY_USER}")


def sa_task_logout():
    SA_TASK_SESSION.delete(f"{ENDPOINT}{app_api_base}/authn?user={CELERY_USER}")
    SA_TASK_SESSION.close()


@celery_app.task(name="sa_task", acks_late=True,
                 # default_retry_delay=10,
                 # autoretry_for=(Exception,),
                 # retry_kwargs={'max_retries': 5},
                 # retry_backoff=True,
                 )
def sa_task(process: str):
    """
    Launch a specific sysadmin task
    :param process:
    :return:
    """

    try:
        sa_task_login()
        print('SA_TASK: ' + process)
        from importlib import import_module
        print(import_module(SA_PROC_PATH + process).run())
        sa_task_logout()
        return None
    except Exception as e:
        log_exception(e)
        return None


@celery_app.task(name="sa_seq_ann_task", acks_late=True,
                 # default_retry_delay=10,
                 # autoretry_for=(Exception,),
                 # retry_kwargs={'max_retries': 5},
                 # retry_backoff=True,
                 )
def sa_seq_ann_task(ann_entries: dict):
    """
    Async import of sequence annotations
    """

    try:
        sa_task_login()
        print('SA_TASK: sequences import_annotations')
        from .sysadmin.sequences.import_annotations import run
        run(ann_entries)
        sa_task_logout()
        return None
    except Exception as e:
        log_exception(e)
        return None
