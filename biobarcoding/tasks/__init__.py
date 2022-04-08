import socket

from celery import Celery

import biobarcoding as base_app_pkg
from .celeryconfig import celery_config


def is_port_open(host="localhost", port=6379):
    a_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    location = (host, port)
    try:
        result_of_check = a_socket.connect_ex(location)
    except:
        result_of_check = 1
    a_socket.close()

    return result_of_check == 0


def get_redis_host():
    if is_port_open("localhost", 6379):
        host = "localhost"
    elif is_port_open("redis", 6379):
        host = "redis"
    else:
        # import redislite
        # redislite.Redis(host="localhost")  # serverconfig={'port': '6379'}
        # if is_port_open("localhost", 6379):
        #     host = "localhost"
        # else:
        raise Exception("No REDIS instance available")

    return host


def initialize_celery(flask_app):
    if "CELERY_BROKER_URL" in flask_app.config:
        broker_url = flask_app.config["CELERY_BROKER_URL"]
    else:
        host = get_redis_host()
        broker_url = f'redis://{host}:6379/'

    if "CELERY_BACKEND_URL" in flask_app.config:
        backend_url = flask_app.config["CELERY_BACKEND_URL"]
    else:
        host = get_redis_host()
        backend_url = f'redis://{host}:6379/'

    # create context tasks in celery
    celery = Celery(
        flask_app.import_name,
        broker=broker_url,
        backend=backend_url
    )
    TaskBase = celery.Task

    class ContextTask(TaskBase):
        abstract = True

        def __call__(self, *args, **kwargs):
            with flask_app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)

    celery.Task = ContextTask

    base_app_pkg.celery = celery


# CELERY APP, used by CELERY TASKS
celery_app = base_app_pkg.celery

if not celery_app:
    host = get_redis_host()
    celery_app = Celery('__main__', broker=f'redis://{host}:6379/', backend=f'redis://{host}:6379/')
    celery_app.autodiscover_tasks(["biobarcoding.tasks.definitions"])
    celery_app.conf.update(celery_config)
