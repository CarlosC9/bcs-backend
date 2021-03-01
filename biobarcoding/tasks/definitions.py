import json
import os

import requests
from billiard.context import Process

from biobarcoding.jobs.ssh_resource import JobExecutorWithSSH
from biobarcoding.tasks import celery_app
from time import sleep, time

from biobarcoding.common.decorators import celery_wf
from biobarcoding.common import check_pid_running
from biobarcoding.jobs import JobExecutorAtResourceFactory
from biobarcoding.db_models.jobs import Job
from biobarcoding.common import ROOT

"""
Celery tasks CANNOT be debugged in Celery!! (the code is run in a separate process; of course they can be debugged "off-line")
"""


# Send messages
# Refresh queues of jobs (at different computing resources)


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


# "job_context" for Celery tasks:
#
# "job_id": ...,
# "endpoint_url": ...
# "resource": {
#     "name": "",
#     "jm_type": "",
#     "jm_location": {
#     },
#     "jm_credentials": {
#     }
# }
# "process": {
#     "name": "",
#     "inputs": {
#     }
# }

@celery_app.task(name="prepare")
@celery_wf(celery_app, wf1, "prepare")
def wf1_prepare_workspace(job_context):
    """
    Prepare workspace for execution of a Job
    :param job_context:
    :return:
    """
    # tmp = json.loads(job_context)

    # Example access RESTful endpoint of "bcs-backend"
    # requests.get(job_context["endpoint_url"]+f"/api/jobs/{job_context['job_id']}")

    # TODO update Job status to "preparing_workspace"
    #  requests.put(job_context["endpoint_url"] + f"/api/jobs/{job_context['job_id']}/status", "preparing_workspace")

    # Get resource manager: ssh, galaxy, other
    # job_executor = JobExecutorAtResourceFactory()
    # job_executor = JobExecutorAtResourceFactory.get(tmp["resource"]["jm_type"], tmp["process"])

    # Launch subprocess if needed
    # if "_pid" not in tmp:
    #     p = Process(target=dummy_func, args=(outfile, 13))
    #     p.start()
    #     tmp["_pid"] = p.pid
    #     append_text(outfile, f"prepare_workspace. PID: {p.pid}")
    #     job_context = json.dumps(tmp)
    #
    # # TODO Task "work"
    # append_text(outfile, "prepare_workspace")
    #
    # # Wait for subprocess to finish (if needed)
    # if check_pid_running(tmp["_pid"]):
    #     append_text(outfile, "retrying task ...")
    #     return 3, job_context  # Return a tuple with first element <seconds to wait>, <context> to repeat the task
    # else:  # Clear "_pid" if subprocess was launched
    #     del tmp["_pid"]
    #     job_context = json.dumps(tmp)
    #     append_text(outfile, "prepare_workspace FINISHED")
    #     return job_context  # Return nothing (None) or <context> (if context changed) to move to the default next task
    # TODO update Job status to "preparing_workspace"
    #  requests.put(job_context["endpoint_url"] + f"/api/jobs/{job_context['job_id']}/status", "preparing_workspace")
    pass


@celery_app.task(name="export")
@celery_wf(celery_app, wf1, "export")
def wf1_export_to_supported_file_formats(job_context: str):
    """
    Prepare input files by exporting data using RESTful services

    :param job_context:
    :return:
    """
    # tmp = json.loads(job_context)

    # Example access RESTful endpoint of "bcs-backend"
    # requests.get(job_context["endpoint_url"]+f"/api/jobs/{job_context['job_id']}")

    # TODO update Job status to "preparing_workspace"
    #  requests.put(job_context["endpoint_url"] + f"/api/jobs/{job_context['job_id']}/status", "preparing_workspace")

    # Get resource manager: ssh, galaxy, other
    # job_executor = JobExecutorAtResourceFactory.get(tmp["resource"]["jm_type"], tmp["process"])

    append_text("export_to_supported_file_formats")
    sleep(2)


@celery_app.task(name="transfer_data")
@celery_wf(celery_app, wf1, "transfer_data")
def wf1_transfer_data_to_resource(job_context: object) -> object:
    """
    Transfer data to the compute resource
    :param job_context:
    :return:
    """
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../tests/ssh_data_test")
    filepackage_list = [dict(transfer_type="upload", local_path=os.path.join(data_dir, "myfile.txt"),
                             remote_path="myfile.txt"),
                        dict(transfer_type="upload", local_path=os.path.join(data_dir, "add_line.py"),
                             remote_path="add_line.py"),
                        dict(transfer_type="upload", local_path=os.path.join(data_dir, "count_lines.py"),
                             remote_path="count_lines.py"),
                        dict(transfer_type="upload", local_path=os.path.join(data_dir, "test.sh"),
                             remote_path="test.sh")]

    tmp = json.loads(job_context)
    if "upload_files" not in tmp["process"]["inputs"]["parameters"]:
        tmp["process"]["inputs"]["parameters"]["upload_files"] = filepackage_list

    job_executor = JobExecutorAtResourceFactory().get(tmp["resource"], tmp["process"]["inputs"]["parameters"])

    transfer_state = tmp.get("transfer_state")
    print(f"Transfer state: {transfer_state}")
    if transfer_state:
        i = transfer_state["idx"]
        pid = transfer_state["pid"]
    else:
        i = 0
        pid = None

    # Ith transfer
    transfer_at_i = tmp["process"]["inputs"]["parameters"]["upload_files"][i] if i < len(filepackage_list) else dict(local_path="", remote_path="")
    local_path = transfer_at_i["local_path"]
    remote_path = transfer_at_i["remote_path"]  # Add workspace base?

    if i == len(filepackage_list):  # Transfer finished
        print("Transfer finished")
        del tmp["transfer_state"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.job_status(pid) == "running":  # Transfer is being executed
        print("Transfer executing")
        sleep(5)
        return None, job_context
    elif job_executor.exists(local_path, remote_path):  # File i has been transferred successfully
        print(f"File {i} transferred: : {local_path} -> {remote_path}. Moving to next")
        i += 1
        tmp["transfer_state"] = dict(idx=i, pid=None)
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i}: {local_path} -> {remote_path}")
        pid = job_executor.upload_file("", local_path, remote_path)
        tmp["transfer_state"] = dict(idx=i, pid=pid)
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
    job_executor = JobExecutorAtResourceFactory()
    print(tmp["process"]["inputs"]["parameters"])
    job_executor = job_executor.get(tmp["resource"], tmp["process"]["inputs"]["parameters"])
    #TODO: ver galaxy -- inputs = tmp["process"]
    inv_id = job_executor.submit(str(tmp['job_id']), tmp["process"]["inputs"]["parameters"])
    tmp['g_id'] = inv_id
    job_context = json.dumps(tmp)
    append_text(f"submit. workspace: {tmp['g_id']}")
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
    job_executor = JobExecutorAtResourceFactory()
    job_executor = job_executor.get(tmp["resource"]["jm_type"], tmp['resource'])
    status = job_executor.job_status(tmp["g_id"])
    if status == 'running':  # pasa siempre or running?
        append_text(f"wait_until_execution_starts: status: {status}")
        return job_context
    if status == 'error':
        append_text(f"wait_until_execution_starts: status: {status}")
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
    job = JobExecutorAtResourceFactory()
    job_executor = job.get(tmp["resource"]["jm_type"], tmp["resource"])
    status = job_executor.job_status(tmp["g_id"])
    if isinstance(status, dict):
        append_text(f"wait_for_execution_end: status: {status}")
        return 'error', job_context
    if status == 'ok':
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
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory()
    job_executor = job_executor.get(tmp["resource"], tmp["process"]["inputs"]["parameters"])
    r = job_executor.get_results(tmp["g_id"])
    tmp['results'] = r
    job_context = json.dumps(tmp)
    append_text(f"transfer_data_from: results: {r}")
    return job_context


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
    job_executor = JobExecutorAtResourceFactory()
    job_executor = job_executor.get(tmp["resource"], tmp["process"]["inputs"]["parameters"])
    job_executor.remove_job_workspace(str(tmp["job_id"]))
    del tmp['g_id']
    job_context = json.dumps(tmp)
    append_text(f"cleanup:")
    return job_context


@celery_app.task(name="success")
@celery_wf(celery_app, wf1, "success")
def wf1_complete_succesfully(job_context: str):
    """
    Just mark the Job as "completed succesfully"

    :param job_context:
    :return:
    """
    append_text("complete_successfully")
    sleep(2)


@celery_app.task(name="error")
@celery_wf(celery_app, wf1, "error")
def wf1_completed_error(job_context: str):
    """
    Mark the Job as "completed with error"

    :param job_context:
    :return:
    """

    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory()
    job_executor = job_executor.get(tmp["resource"], tmp["process"]["inputs"]["parameters"])
    status = job_executor.job_status(job_executor["g_id"])
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
    job_executor = JobExecutorAtResourceFactory()
    job_executor = job_executor.get(tmp["resource"], tmp["process"]["inputs"]["parameters"])
    job_executor.cancel_job(tmp['g_id'])
    job_executor.remove_job_workspace(str(tmp['job_id']))
    append_text(f"cancelled")
    return job_context
