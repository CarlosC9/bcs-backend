import json
import os

import requests
from billiard.context import Process

from biobarcoding.rest.file_manager import FilesAPI
from biobarcoding.tasks import celery_app
from time import sleep, time

from biobarcoding.common.decorators import celery_wf
from biobarcoding.common import check_pid_running
from biobarcoding.jobs import JobExecutorAtResourceFactory
from biobarcoding.db_models.jobs import Job
from biobarcoding.common import ROOT
from biobarcoding.rest import bcs_api_base


"""
Celery tasks CANNOT be debugged in Celery!! (the code is run in a separate process; of course they can be debugged "off-line")
"""

MAX_ERRORS = 3

# Send messages
# Refresh queues of jobs (at different computing resources)

def change_status(tmp, status: str):
    endpoint_url, job_id = tmp["endpoint_url"], tmp["job_id"]
    url = f"{endpoint_url}/api/jobs/{job_id}"
    if tmp["status"] != status:
        # todo pasar el cookie en el job context
        cmd = ["curl", "-c", "cookies.txt", "-X", "PUT", "http://localhost:5000/api/authn?user=test_user", "&&", "curl",
               "-b", "cookies.txt", "-X", "PUT", "-H", "Content-Type: application/json", "-d",
               '{"status":"preparing_workspace"}', url]
        import subprocess
        print(subprocess.run(cmd))
        tmp["status"] = status
        return tmp
    else:
        return tmp

@celery_app.task
def add(x, y):
    print(f"NORMAL CELERY TASK: {x}+{y}={x + y}")
    return x + y


@celery_app.task
def periodic_sum(x, y):
    print(f"PERIODIC CELERY TASK: {x}+{y}={x + y}")
    return x + y


# ----------------------------------------------------------------

def append_text(s: str):
    file = ROOT + "/tests/data_test/celery_log.txt"
    with open(file, "a+") as f:
        f.write(f"{s}\n")


# To test:
# - Run "bcs-backend" in the IDE
# - Open three shells in the "bcs-backend" directory:
#   - ./start_local_dev_services.sh
#   - python3 biobarcoding/rest/main.py
#   - celery -A biobarcoding.tasks.definitions flower [OPTIONAL: to see tasks, http://localhost:5555]
#   - curl -i -XPOST http://localhost:5000/api/jobs/ --data-urlencode "{}" --> (should return immediately)
#     - (to see messages generated by tasks of this workflow) tail -f /home/rnebot/Downloads/borrame/log.txt
#   - curl -i -XPOST http://localhost:5000/api/jobs/ -H "Content-Type: application/json" -d @"/home/daniel/Documentos/GIT/bcs-backend/tests/request_transfer.json"
wf1 = {
    "prepare": {"": "export"},
    "export": {"": "transfer_data"},
    "transfer_data": {"": "submit"},
    "submit": {"": "wait_until_execution_starts"},
    "wait_until_execution_starts": {
        "": "wait_for_execution_end",
        "error": "error"

    },
    "wait_for_execution_end": {
        "": "transfer_data_from",
        "error": "error"

    },
    "transfer_data_from": {"": "import"},
    "import": {"": "cleanup"},
    "cleanup": {"": "success"},
    # "success"
    # "error"
    # "cancel"
}


def dummy_func(file, secs):
    print(f"Dummy {file}, {secs}")
    start = time()
    endc = start
    while (endc - start) < secs:
        append_text(f"elapsed {endc - start} of {secs}")
        sleep(6)
        endc = time()
    os.exit(0)



job_context = {"endpoint_url": "http//:localhost:5000/",
               "process":
                   {"inputs":
                        {"parameters": "...."},
                         "data": [{
                                "bo_type": "collection",
                                "ids": ["jkdjdjdjjdjd"],
                                "type": "fasta"},
                                  {
                                    "bo_type": "path",
                                      "ids":["...."],
                                      "type":"fasta"
                                  },
                                 {"so_type": "seq",
                                  "bo": ["jkdjdjdjjdjd", "jdjdjdjfjdjdj"],
                                "type": "fasta"}
                         ]
                    },
                    "name": {"........"},
               "outputs":[{"....."}],
               "status": "created",
               "resource": {"......."},
               "job_id": 60}






@celery_app.task(name="prepare")
@celery_wf(celery_app, wf1, "prepare")
def wf1_prepare_workspace(job_context):
    """
    Prepare workspace for execution of a Job
    :param job_context:
    :return:
    """
    tmp = json.loads(job_context)
    tmp = change_status(tmp,"preparing_workspace")

    # Example access RESTful endpoint of "bcs-backend"
    # requests.get(job_context["endpoint_url"]+f"/api/jobs/{job_context['job_id']}")

    # TODO update Job status to "preparing_workspace"
    #  requests.put(job_context["endpoint_url"] + f"/api/jobs/{job_context['job_id']}/status", "preparing_workspace")

    # Get resource manager: ssh, galaxy, other
    # job_executor = JobExecutorAtResourceFactory()
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    job_executor.create_job_workspace(str(tmp["job_id"]))


@celery_app.task(name="export")
@celery_wf(celery_app, wf1, "export")
def wf1_export_to_supported_file_formats(job_context: str):
    """
    Prepare input files by exporting data using RESTful services

    :param job_context:
    :return:
    """
    # tmp = json.loads(job_context)
    #tmp = change_status(tmp, "preparing_workspace")

    # Example access RESTful endpoint of "bcs-backend"
    # requests.get(job_context["endpoint_url"]+f"/api/jobs/{job_context['job_id']}")

    # TODO update Job status to "preparing_workspace"
    #  requests.put(job_context["endpoint_url"] + f"/api/jobs/{job_context['job_id']}/status", "preparing_workspace")

    # Get resource manager: ssh, galaxy, other
    # job_executor = JobExecutorAtResourceFactory.get(tmp["resource"]["jm_type"], tmp["process"])

    append_text("export_to_supported_file_formats")
    sleep(2)

job_context = {"endpoint_url": "http//:localhost:5000/",
               "process":
                   {"inputs":
                        {"parameters": {"ClustalW": {"darna": "PROTEIN"}},
                         "data": [
                             {"path": "...."},
                             { "path": "path"},
                             {"type": "fasta"}
                         ]},
                    "name": "MSA ClustalW"},
               "status": "created",
               "resource": {".........."},
               "job_id": 60}



@celery_app.task(name="transfer_data")
@celery_wf(celery_app, wf1, "transfer_data")
def wf1_transfer_data_to_resource(job_context: object) -> object:
    """
    Transfer data to the compute resource
    :param job_context:
    :return:
    """
    tmp = json.loads(job_context)
    tmp = change_status(tmp, "transfer_data_to_resource")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    transfer_state = tmp.get("transfer_state")
    print(f"Transfer state: {transfer_state}")
    if transfer_state:  # ya ha empezado la transferencia
        i = transfer_state["idx"]
        pid = transfer_state["pid"]  # job_id
        n_errors = transfer_state["n_errors"]
    else:
        i = 0
        pid = None
        n_errors = 0
        tmp["transfer_state"] = dict(idx=i, pid=None, n_errors=n_errors)

    # Ith transfer
    files_list = tmp["process"]["inputs"]["data"]
    transfer_at_i = files_list[i] if i < len(files_list) else dict(
        path="")

    '''The transfer_at_i keys are unique for each job type. The only field that is mandatory
    is the path that refers to the local path of the file to be transferred'''
    local_path = transfer_at_i.get("path")

    if i < len(files_list) and not os.path.isfile(local_path):
        tmp["transfer_state"] = f"{local_path} doesn't exist"
        job_context = json.dumps(tmp)
        return "error", job_context
    elif i == len(files_list):  # Transfer finished
        print("Transfer finished")
        del tmp["transfer_state"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.job_status(pid) == "running":  # Transfer is being executed
        print("Transfer executing")
        sleep(5)
        return None, job_context
    elif job_executor.exists(tmp):  # File i has been transferred successfully
        print(f"File {i} transferred: {local_path} . Moving to next")
        i += 1
        tmp["transfer_state"] = dict(idx=i, pid=None, n_errors=n_errors)
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i}: {local_path}")
        pid = job_executor.upload_file(tmp)
        # TODO yo puedo tener un error aquí
        tmp["transfer_state"] = dict(idx=i, pid=pid, n_errors=n_errors)
        print(tmp['transfer_state'])
        job_context = json.dumps(tmp)
        return None, job_context


    
@celery_app.task(name="submit")
@celery_wf(celery_app, wf1, "submit")
def wf1_submit(job_context: str):
    """
    Submit job to compute resource
    :param job_context:
    :return:
    """
    tmp = json.loads(job_context)
    tmp = change_status(tmp, "submit")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    inputs = tmp["process"]
    pid = job_executor.submit(tmp["process"])
    tmp['pid'] = pid
    job_context = json.dumps(tmp)
    append_text(f"submit. workspace: {tmp['pid']}")
    return job_context


@celery_app.task(name="wait_until_execution_starts")
@celery_wf(celery_app, wf1, "wait_until_execution_starts")
def wf1_wait_until_execution_starts(job_context: str):
    """
    Wait for the job to start executing at the
    :param job_context:
    :return:
    """
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    status = job_executor.job_status(tmp["pid"])
    if status == 'running' or 'ok':  # pasa siempre or running?
        append_text(f"wait_until_execution_starts: status: {status}")
        return job_context
    if status == 'error':
        append_text(f"wait_until_execution_starts: status: {status}")
        # pasar error en job_context?
        return 'error', job_context
    else:
        append_text(f"wait_until_execution_starts: status: {status}")
        return 3, job_context


@celery_app.task(name="wait_for_execution_end")
@celery_wf(celery_app, wf1, "wait_for_execution_end")
def wf1_wait_for_execution_end(job_context: str):
    """
    Wait for the job to finish execution, knowing it is running
    When finished, it can end successfully or with an error.

    :param job_context:
    :return:
    """
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    status = job_executor.job_status(tmp["pid"])
    if isinstance(status, dict):
        append_text(f"wait_for_execution_end: status: {status}")
        return 'error', job_context
    elif status == 'ok':
        append_text(f"wait_for_execution_end: status: {status}")
        return job_context
    else:
        append_text(f"wait_for_execution_end: status: {status}")
        return 3, job_context


@celery_app.task(name="transfer_data_from")
@celery_wf(celery_app, wf1, "transfer_data_from")
def wf1_transfer_data_from_resource(job_context: str):
    """
    Once the computation ends, transfer results from the resource to local
    :param job_context:
    :return:
    """

    """
        Transfer data to the compute resource
        :param job_context:
        :return:
        """
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    transfer_state = tmp.get("transfer_state")
    print(f"Transfer state: {transfer_state}")
    if transfer_state:  # ya ha empezado la transferencia
        i = transfer_state["idx"]
        pid = transfer_state["pid"]  # job_id
    else:
        i = 0
        pid = None
        os.mkdir(f"base_path_to_results/{tmp['job_id']}/")
        tmp["transfer_state"] = dict(idx=i, pid=None)

    # Ith transfer
    files_list = tmp["process"]["inputs"]["results"] # en mi caso esta es una lista de ids que tengo que pedir a mi workspace
    transfer_at_i = files_list[i] if i < len(files_list) else dict(
        local_path="", remote_path="")

    '''The transfer_at_i keys are unique for each job type. The only field that is mandatory
    is the path that refers to the local path of the file to be transferred'''
    local_path = os.path.join(f"base_path_to_results/{tmp['job_id']}/", transfer_at_i.get("path"))

    # miss transferencias no tienen estados
    if i == len(files_list):  # Transfer finished
        print("Transfer finished")
        del tmp["transfer_state"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.job_status(pid) == "running":  # Transfer is being executed
        print("Transfer executing")
        sleep(5)
        return None, job_context
    elif os.path.isfile(local_path):  # File i has been transferred successfully
        print(f"File {i} transferred: {local_path} -> Moving to next")
        FilesAPI.put(local_path)
        i += 1
        tmp["transfer_state"] = dict(idx=i, pid=None)
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i}: {local_path}")
        pid = job_executor.download_file(tmp)
        # TODO yo puedo tener un error aquí
        tmp["transfer_state"] = dict(idx=i, pid=pid)
        print(tmp['transfer_state'])
        job_context = json.dumps(tmp)
        return None, job_context


@celery_app.task(name="import")
@celery_wf(celery_app, wf1, "import")
def wf1_import_into_database(job_context: str):
    """
    Import resulting files into the database

    :param job_context:
    :return:
    """
    append_text("import_into_database")
    sleep(2)


@celery_app.task(name="cleanup")
@celery_wf(celery_app, wf1, "cleanup")
def wf1_cleanup_workspace(job_context: str):
    """
    Delete workspace at the remote resource

    :param job_context:
    :return:
    """

    tmp = json.loads(job_context)
    tmp = change_status(tmp, "cleanup")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    job_executor.remove_job_workspace(str(tmp["job_id"]))
    del tmp['pid']
    job_context = json.dumps(tmp)
    append_text(f"cleanup:")


@celery_app.task(name="success")
@celery_wf(celery_app, wf1, "success")
def wf1_complete_succesfully(job_context: str):
    """
    Just mark the Job as "completed succesfully"

    :param job_context:
    :return:
    """
    # TODO CHANGE STATE
    append_text("complete_successfully")
    return job_context


@celery_app.task(name="error")
@celery_wf(celery_app, wf1, "error")
def wf1_completed_error(job_context: str):
    """
    Mark the Job as "completed with error"

    :param job_context:
    :return:
    """

    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    status = job_executor.job_status(job_executor["pid"])
    tmp['error'] = status
    job_executor.remove_job_workspace(str(tmp['job_id']))
    job_context = json.dumps(tmp)
    append_text(f"error: {tmp['error']}")
    return job_context


@celery_app.task(name="cancel")
@celery_wf(celery_app, wf1, "cancel")
def wf1_cancelled(job_context: str):
    """
    Mark the Job as "cancelled"

    :param job_context:
    :return:
    """
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    job_executor.cancel_job(tmp['pid'])
    job_executor.remove_job_workspace(tmp)
    append_text(f"cancelled")
    return job_context
