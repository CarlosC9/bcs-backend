import json
import os

import requests
from billiard.context import Process

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

def append_text( s: str):
    file = ROOT + "/tests/data_test/celery_log.txt"
    with open(file, "a+") as f:
        f.write(f"{s}\n")


# To test:
# - Run "bcs-backend" in the IDE
# - Open three shells in the "bcs-backend" directory:
#   - ./start_local_dev_services.sh
#   - celery -A biobarcoding.tasks.definitions flower [OPTIONAL: to see tasks, http://localhost:5555]
#   - curl -i -XPOST http://localhost:5000/api/jobs/ --data-urlencode "{}" --> (should return immediately)
#     - (to see messages generated by tasks of this workflow) tail -f <project folder>/tests/data_test/celery_log.txt

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
    tmp = json.loads(job_context)

    # Example access RESTful endpoint of "bcs-backend"
    # requests.get(job_context["endpoint_url"]+f"/api/jobs/{job_context['job_id']}")

    # TODO update Job status to "preparing_workspace"
    #  requests.put(job_context["endpoint_url"] + f"/api/jobs/{job_context['job_id']}/status", "preparing_workspace")

    # Get resource manager: ssh, galaxy, other
    job_executor = JobExecutorAtResourceFactory.get(tmp["resource"]["jm_type"], tmp["process"])
    wid = job_executor.create_job_workspace(tmp['job_id'])
    tmp['process']['w_id'] = wid
    job_context = json.dumps(tmp)
    # se puede dar que ese workspace ya exista o no funcione galaxy

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
    return job_context


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
def wf1_transfer_data_to_resource(job_context: str):
    """
    Transfer data to the compute resource
    :param job_context:
    :return:
    """
    append_text("transfer_data_to_resource")
    sleep(2)

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
    job_executor = job_executor.get(tmp["resource"]["jm_type"], tmp['resource'])
    inputs = tmp["process"]
    inv_id = job_executor.submit(str(tmp['job_id']), inputs)
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
    job_executor = job_executor.get(tmp["resource"]["jm_type"], tmp["resource"])
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
    job_executor = job_executor.get(tmp["resource"]["jm_type"], tmp["resource"])
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
    job_executor = job_executor.get(tmp["resource"]["jm_type"], tmp["resource"])
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
    job_executor = job_executor.get(tmp["resource"]["jm_type"], tmp["resource"])
    job_executor.cancel_job(tmp['g_id'])
    job_executor.remove_job_workspace(str(tmp['job_id']))
    append_text(f"cancelled")
    return job_context