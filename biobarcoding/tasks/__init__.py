import socket

from celery import Celery
from sqlalchemy import and_

import biobarcoding
from .celeryconfig import celery_config
from ..db_models import DBSession
from ..db_models.sysadmin import Identity,Authenticator,IdentityAuthenticator
from ..rest import bcs_api_base


def is_port_open(host="localhost", port=6379):
    a_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    location = (host, port)
    try:
        result_of_check = a_socket.connect_ex(location)
    except:
        result_of_check = 1
    a_socket.close()

    return result_of_check == 0

def create_celery_user():
    session = DBSession()
    test_user_id = "celery_user"
    authenticator_id = "75f3373b-ee0e-4e13-85da-2045010f3939"
    iden = session.query(Identity).filter(Identity.name == test_user_id).first()
    authentication = session.query(Authenticator).filter(Authenticator.uuid == authenticator_id).first()
    iden_authentication = session.query(IdentityAuthenticator). \
        filter(
        and_(IdentityAuthenticator.identity == iden, IdentityAuthenticator.authenticator == authentication)).first()
    if not iden_authentication:
        iden_authentication = IdentityAuthenticator()
        iden_authentication.identity = iden
        iden_authentication.authenticator = authentication
        iden_authentication.name = "celery_user"
        iden_authentication.email = "celery@celery.org"
        session.add(iden_authentication)
        #TODO rol: only put job
    session.commit()



def change_status(status,job_context):
    import json
    import requests
    tmp = json.loads(job_context)
    job_id = tmp["job_id"]
    if tmp["status"] != status:
        url = tmp["endpoint"]
        s = requests.Session()
        # todo permisos propios para celery
        response = s.put(f"{url}{bcs_api_base}/authn?user=test_user")
        print(f"status code: {response.status_code}")
        response = s.put(f'{url}{bcs_api_base}/jobs/{job_id}', json={'status': status})
        print(f"status code: {response.status_code}")
        return response.status_code
    else:
        return None


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

    biobarcoding.celery = celery
#    create_celery_user()


# CELERY APP, used by CELERY TASKS
celery_app = biobarcoding.celery

if not celery_app:
    host = get_redis_host()
    celery_app = Celery('__main__', broker=f'redis://{host}:6379/', backend=f'redis://{host}:6379/')
    celery_app.autodiscover_tasks(["biobarcoding.tasks.definitions"])
    celery_app.conf.update(celery_config)
