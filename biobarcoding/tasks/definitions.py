from . import celery_app
from time import sleep

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

"""
 "created",
 "preparing_workspace",
 "exporting_to_supported_file_formats",
 "transferring_data_to_resource",
 "submitting",
 "waiting_for_execution",
 "executing",
   "transferring_data_from_resource",
     "importing_into_database",
       "cleaning_up_workspace",
         "completed_successfully",
   "completed_error",
   "cancelled",

job_context:
{
  "job_id": "..",
  "endpoint_url": "<base to call>",
  "resource": {
    "name": "",
    "proc_arch": "",
    "operating_system": "",
    "jm_type": "",
    "jm_location": {
    },
    "jm_credentials": {
        
    }
  }
  "workflow": {
    "name": "",
    "id": "",
    "inputs": {
    }
  }
"""

outfile = "/home/rnebot/Downloads/borrame/log.txt"


def append_text(file: str, s: str):
    with open(file, "a+") as f:
        f.write(f"{s}\n")


@celery_app.task
def wf1_prepare_workspace(job_context: str):
    """
    Prepare workspace for execution of a Job
    :param job_context:
    :return:
    """
    # TODO Obtain RESTful endpoint.
    #  Read current status and, if it is "created",
    #  update Job status to "preparing_workspace"

    # TODO Read resource type: ssh, galaxy, other
    # TODO Create resource manager instance and call
    append_text(outfile, "prepare_workspace")
    sleep(2)


@celery_app.task
def wf1_export_to_supported_file_formats(job_context: str):
    """
    Prepare input files by exporting data using RESTful services

    :param job_context:
    :return:
    """
    # TODO
    append_text(outfile, "export_to_supported_file_formats")
    sleep(2)


@celery_app.task
def wf1_transfer_data_to_resource(job_context: str):
    """
    Transfer data to the compute resource
    :param job_context:
    :return:
    """
    append_text(outfile, "transfer_data_to_resource")
    sleep(2)


@celery_app.task
def wf1_submit(job_context: str):
    """
    Submit job to compute resource
    :param job_context:
    :return:
    """
    append_text(outfile, "submit")
    sleep(2)


@celery_app.task
def wf1_wait_until_execution_starts(job_context: str):
    """
    Wait for the job to start executing at the
    :param job_context:
    :return:
    """
    append_text(outfile, "wait_until_execution_starts")
    sleep(2)


@celery_app.task
def wf1_wait_for_execution_finish(job_context: str):
    """
    Wait for the job to start executing at the
    :param job_context:
    :return:
    """
    append_text(outfile, "wait_for_execution_to_finish")
    sleep(2)


@celery_app.task
def wf1_transfer_data_from_resource(job_context: str):
    """
    Once the computation ends, transfer results from the resource to local
    :param job_context:
    :return:
    """
    append_text(outfile, "transfer_data_from_resource")
    sleep(2)


@celery_app.task
def wf1_import_into_database(job_context: str):
    """
    Import resulting files into the database

    :param job_context:
    :return:
    """
    append_text(outfile, "import_into_database")
    sleep(2)


@celery_app.task
def wf1_cleanup_workspace(job_context: str):
    """
    Delete workspace at the remote resource

    :param job_context:
    :return:
    """
    append_text(outfile, "cleanup_workspace")
    sleep(2)


@celery_app.task
def wf1_complete_succesfully(job_context: str):
    """
    Just mark the Job as "completed succesfully"

    :param job_context:
    :return:
    """
    append_text(outfile, "complete_successfully")
    sleep(2)


@celery_app.task
def wf1_completed_error(job_context: str):
    """
    Mark the Job as "completed with error"

    :param job_context:
    :return:
    """
    append_text(outfile, "completed_error")
    sleep(2)


@celery_app.task
def wf1_cancelled(job_context: str):
    """
    Mark the Job as "cancelled"

    :param job_context:
    :return:
    """
    append_text(outfile, "cancelled")
    sleep(2)
