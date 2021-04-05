import json
import os
import requests
import subprocess

from biobarcoding.rest import bcs_api_base
from biobarcoding.tasks import celery_app
from time import sleep

from biobarcoding.common.decorators import celery_wf
from biobarcoding.jobs import JobExecutorAtResourceFactory
from biobarcoding.common import ROOT

"""
Celery tasks CANNOT be debugged in Celery!! (the code is run in a separate process; of course they can be debugged "off-line")
"""

MAX_ERRORS = 3
TMP_PATH = "/tmp"


# Send messages
# Refresh queues of jobs (at different computing resources)

def change_status(tmp, status: str):
    endpoint_url = os.getenv("ENDPOINT_URL")
    cookies_file_path = os.getenv("COOKIES_FILE_PATH")
    print(f"Env. variable ENDPOINT_URL: {endpoint_url}")
    print(f"Env. variable COOKIES_FILE_PATH: {cookies_file_path}")
    job_id = tmp["job_id"]
    # TODO comprobar existencia de
    url = f"{endpoint_url}{bcs_api_base}/jobs/{job_id}"
    if tmp["status"] != status:
        api_login()
        #TODO: este comando hace que se borre la sesión. Hay que revisar que esté funcionando bien.
        cmd = ["curl", "--cookie-jar", cookies_file_path, "-c", cookies_file_path, "-X", "PUT", "-H",
               "Content-Type: application/json", "-d", '{"status":"preparing_workspace"}', url]
        print(subprocess.run(cmd))
        print_file(cookies_file_path)
        tmp["status"] = status
        return tmp
    else:
        return tmp

def print_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            print(line)

def api_login():
    endpoint_url = os.getenv("ENDPOINT_URL")
    cookies_file_path = os.getenv("COOKIES_FILE_PATH")
    url = f"{endpoint_url}{bcs_api_base}/authn?user=test_user"
    cmd = ["curl", "--cookie-jar", cookies_file_path, "-X", "PUT", url]
    subprocess.run(cmd)
    print_file(cookies_file_path)



def check_file_is_stored_in_backend(filename, job_id):
    endpoint_url = os.getenv("ENDPOINT_URL")
    cookies_file_path = os.getenv("COOKIES_FILE_PATH")
    print(f"Env. variable ENDPOINT_URL: {endpoint_url}")
    print(f"Env. variable COOKIES_FILE_PATH: {cookies_file_path}")
    api_login()
    cmd = ["curl", "--cookie-jar", cookies_file_path, "--cookie",
           cookies_file_path, f"{endpoint_url}{bcs_api_base}/files/jobs/{job_id}/{filename}"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    print_file(cookies_file_path)
    process_return_dict = json.loads(proc.stdout)
    if process_return_dict.get("content"):
        return True
    else:
        return False


#TODO: Hay que prepararlo para las colecciones y los ficheros que suba el usuario de manera
# que se puedan concatenar si se refieren al mismo fichero.
# Mirar: https://stackoverflow.com/questions/40359012/how-to-append-a-file-with-the-existing-one-using-curl

def export(file_dict, results_dir) -> object:
    """
    Execute remote client script
    @param url: url of the get
    @param path: local path of the file gotten with curl
    @return: pid: PID of the executed script process
    """
    tmp_path = os.path.join(results_dir, file_dict["file"])
    extension = file_dict['type']
    bos_type = file_dict['bo_type']
    selection_dict = {"filter": [{"feature_id": {"op": "in", "unary": file_dict['selection']}}]}
    selection_json = json.dumps(selection_dict)
    endpoint = os.getenv("ENDPOINT_URL")
    cookies_file_path = os.getenv("COOKIES_FILE_PATH")
    api_login()
    url = f"{endpoint}{bcs_api_base}/bos/{bos_type}.{extension}"
    curl = f"curl --cookie-jar {cookies_file_path} --cookie {cookies_file_path} -X GET -d \'{selection_json}\' -H \'Content-Type:application/json\' {url} -o {tmp_path}"
    cmd = f"(nohup bash -c &quot;{curl}&quot; >/tmp/mtest </dev/null 2>/tmp/mtest.err & echo $!; wait $!; echo $?  >> {results_dir}/$!.exit_status)"
    print(cmd)
    popen_pipe = os.popen(cmd)
    print_file(cookies_file_path)
    pid = popen_pipe.readline().rstrip()
    print(f"PID: {pid}")
    return pid


def is_bos_file(input_file):
    if os.path.exists(input_file):  # and check that it is not an error file:
        return True
    else:
        return False


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
#   - curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt -XPOST http://localhost:5000/api/jobs/ -H "Content-Type: application/json" -d @"/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/galaxy_request.json"
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
    "transfer_data_from": {"": "store_result_in_backend"},
    "store_result_in_backend": {"": "cleanup"},
    "cleanup": {"": "success"},
    # "success"
    # "error"
    # "cancel"
}


@celery_app.task(name="prepare")
@celery_wf(celery_app, wf1, "prepare")
def wf1_prepare_workspace(job_context):
    """
    Prepare workspace for execution of a Job
    :param job_context:
    :return:
    """
    '''cookies_file_path = os.getenv("COOKIES_FILE_PATH")
    endpoint_url = os.getenv("ENDPOINT_URL")
    cmd = ["curl", "--cookie-jar", cookies_file_path, f"'{endpoint_url}{bcs_api_base}/authn?user=test_user'"]
    subprocess.run(cmd)'''
    tmp = json.loads(job_context)
    tmp = change_status(tmp, "preparing_workspace")

    job_executor = JobExecutorAtResourceFactory().get(tmp)
    job_executor.create_job_workspace(str(tmp["job_id"]))
    job_context = json.dumps(tmp)
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
    # tmp = change_status(tmp, "preparing_workspace")

    # Example access RESTful endpoint of "bcs-backend"
    # request files from bos ID: 1,2 and save it in local folder
    # curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt -X GET -d '{"filter": [{"feature_id": {"op": "in",
    # "unary": [1, 2]}}]}' -H "Content-Type:application/json"  http://localhost:5000/api/bos/sequences.fasta -o $TEST_FILES_PATH/test.fasta
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    state_dict = tmp.get("state_dict")

    print(f"Export state: {state_dict}")

    if state_dict:  # ya ha empezado la transferencia
        i = state_dict["idx"]
        n_errors = state_dict["n_errors"]
        results_dir = state_dict["results_dir"]
    else:
        i = 0
        n_errors = 0
        results_dir = os.path.join(job_executor.LOCAL_WORKSPACE, str(tmp['job_id']))
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors, state="export", results_dir=results_dir)

    # Files to transfer
    files = tmp['process']['inputs']['data']
    file_dict = files[i] if i < len(files) else {"file": "", "remote_name": "", "type": ""}

    if i < len(files) and n_errors >= MAX_ERRORS:
        print(f"Export error: File {i + 1} {file_dict['file']}")
        job_context = json.dumps(tmp)
        return "error", job_context
    if i == len(files):  # Transfer finished
        print("Export finished")
        del tmp["state_dict"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.job_status(tmp) == "running":  # Transfer is being executed
        print("Export executing")
        return 1, job_context
    elif job_executor.job_status(tmp) == "":
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors + 1, state="export", results_dir=results_dir)
        job_context = json.dumps(tmp)
        return None, job_context
    elif is_bos_file(os.path.join(results_dir, file_dict["file"])):
        print(f"File {file_dict['file']} transferred -> Moving to next")
        tmp['process']['inputs']['data'][i]['file'] = os.path.join(results_dir, file_dict["file"])
        del tmp['process']['inputs']['data'][i]['selection']
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i + 1, n_errors=n_errors, state="export", results_dir=results_dir)
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin export {i + 1}: {file_dict['file']}")
        pid = export(file_dict, results_dir)
        tmp["pid"] = pid
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors, state="export", results_dir=results_dir)
        job_context = json.dumps(tmp)
        return None, job_context


# job_context = {"endpoint_url": "http//:localhost:5000/",
#                "process":
#                    {"inputs":
#                         {"parameters": {"ClustalW": {"darna": "PROTEIN"}},
#                          "data": [
#                              {"path": "...."},
#                              { "path": "path"},
#                              {"type": "fasta"}
#                          ]},
#                     "name": "MSA ClustalW"},
#                "status": "created",
#                "resource": {".........."},
#                "job_id": 60}


@celery_app.task(name="transfer_data")
@celery_wf(celery_app, wf1, "transfer_data")
def wf1_transfer_data_to_resource(job_context: str) -> object:
    """
    Transfer data to the compute resource
    :param job_context:
    :return:
    """
    tmp = json.loads(job_context)
    tmp = change_status(tmp, "transfer_data_to_resource")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    print(tmp)
    state_dict = tmp.get("state_dict")
    print(f"Transfer state: {state_dict}")
    if state_dict:  # ya ha empezado la transferencia
        i = state_dict["idx"]
        n_errors = state_dict["n_errors"]
    else:
        i = 0
        n_errors = 0
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors, state="upload")

    # Files to transfer
    files = job_executor.get_upload_files_list(tmp)
    file_dict = files[i] if i < len(files) else {"file": "", "remote_name": "", "type": ""}

    if i < len(files) and n_errors >= MAX_ERRORS:
        print(f"File {file_dict['file']} doesn't exist.")
        job_context = json.dumps(tmp)
        return "error", job_context
    elif i == len(files):  # Transfer finished
        print("Transfer finished")
        del tmp["state_dict"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.job_status(tmp) == "running":  # Transfer is being executed
        print("Transfer executing")
        return 1, job_context
    elif job_executor.job_status(tmp) == "":  # Transfer error
        print(f"Transfer error: File {i + 1}")
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors + 1, state="upload")
        job_context = json.dumps(tmp)
        return None, job_context
    elif job_executor.exists(tmp):  # File i has been transferred successfully
        print(f"File {file_dict['file']} transferred -> Moving to next")
        i += 1
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors, state="upload")
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i + 1}: {file_dict['file']}")
        pid = job_executor.upload_file(tmp)
        tmp["pid"] = pid
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors, state="upload")
        print(tmp['state_dict'])
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
    status = job_executor.job_status(tmp)
    append_text(f"wait_until_execution_starts: status: {status}")
    if isinstance(status, dict):
        tmp['process']['error'] = status
        job_context = json.dumps(tmp)
        return 'error', job_context
    elif status == 'running' or status == 'ok':
        return job_context
    else:
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
    status = job_executor.job_status(tmp)
    append_text(f"wait_for_execution_end: status: {status}")
    if isinstance(status, dict):
        tmp['process']['error'] = status
        job_context = json.dumps(tmp)
        return 'error', job_context
    elif status == 'ok':
        return job_context
    else:
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
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    state_dict = tmp.get("state_dict")

    print(f"Transfer state: {state_dict}")
    if state_dict:  # ya ha empezado la transferencia
        i = state_dict["idx"]
        n_errors = state_dict["n_errors"]
        results_dir = state_dict["results_dir"]
    else:
        i = 0
        n_errors = 0
        results_dir = os.path.join(job_executor.LOCAL_WORKSPACE, str(tmp['job_id']))
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors, state="download", results_dir=results_dir)

    # Files to transfer
    files = job_executor.get_download_files_list(tmp)
    file_dict = files[i] if i < len(files) else {"file": "", "remote_name": "", "type": ""}

    if i < len(files) and n_errors >= MAX_ERRORS:
        print(f"Transfer error: File {i + 1}")
        job_context = json.dumps(tmp)
        return "error", job_context
    if i == len(files):  # Transfer finished
        print("Transfer finished")
        del tmp["state_dict"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.job_status(tmp) == "running":  # Transfer is being executed
        print("Transfer executing")
        return 1, job_context
    elif job_executor.job_status(tmp) == "":
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors + 1, state="download", results_dir=results_dir)
        job_context = json.dumps(tmp)
        return None, job_context
    elif job_executor.exists(tmp):  # File i has been transferred successfully
        print(f"File {file_dict['file']} transferred -> Moving to next")
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i + 1, n_errors=n_errors, state="download", results_dir=results_dir)
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i + 1}: {file_dict['file']}")
        pid = job_executor.download_file(tmp)
        tmp["pid"] = pid
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors, state="download", results_dir=results_dir)
        job_context = json.dumps(tmp)
        return None, job_context


@celery_app.task(name="store_result_in_backend")
@celery_wf(celery_app, wf1, "store_result_in_backend")
def wf1_store_result_in_backend(job_context: str):
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    state_dict = tmp.get("state_dict")

    print(f"Store result in backend state: {state_dict}")
    if state_dict:  # ya ha empezado la transferencia
        i = state_dict["idx"]
        n_errors = state_dict["n_errors"]
        results_dir = state_dict["results_dir"]
    else:
        i = 0
        n_errors = 0
        results_dir = os.path.join(job_executor.LOCAL_WORKSPACE, str(tmp['job_id']))
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors, results_dir=results_dir, state="store")

    # Files to transfer
    files = job_executor.get_download_files_list(tmp)
    file_dict = files[i] if i < len(files) else {"file": "", "remote_name": "", "type": ""}

    if i < len(files) and n_errors >= MAX_ERRORS:
        print(f"Store result in backend error: File {file_dict['file']}")
        job_context = json.dumps(tmp)
        return "error", job_context
    elif i == len(files):  # Transfer finished
        print("Store result in backend finished")
        del tmp["state_dict"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.job_status(tmp) == "running":  # Transfer is being executed
        print("Store result in backend executing")
        return 1, job_context
    elif job_executor.job_status(tmp) == "":
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors + 1, results_dir=results_dir, state="store")
        job_context = json.dumps(tmp)
        return None, job_context
    elif check_file_is_stored_in_backend(file_dict["file"], tmp["job_id"]):  # File i has been transferred successfully
        print(f"File {file_dict['file']} stored -> Moving to next")
        i += 1
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_errors=n_errors, results_dir=results_dir, state="store")
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Beginning to store result in backend of file {file_dict['file']}")
        cookies_file_path = os.getenv("COOKIES_FILE_PATH")
        endpoint_url = os.getenv("ENDPOINT_URL")
        local_path = os.path.join(job_executor.LOCAL_WORKSPACE, file_dict["file"])
        api_login()
        curl_cmd = f"curl -s --cookie-jar {cookies_file_path} --cookie {cookies_file_path} -H \"Content-Type: application/x-fasta\" -XPUT --data-binary @\"{local_path}\" \"{endpoint_url}{bcs_api_base}/files/jobs/{tmp['job_id']}/{file_dict['file']}.content\""
        # cmd = f"curl -s -c {cookies_file_path} -X PUT {endpoint_url}{bcs_api_base}/authn?user=test_user > /dev/null && (nohup bash -c \'{curl_cmd} \' >/tmp/mtest </dev/null 2>/tmp/mtest.err & echo $!; wait $!; echo $? >> /tmp/{tmp['job_id']}/$!.exit_status)"
        cmd = f"(nohup bash -c \'{curl_cmd} \' >/tmp/mtest </dev/null 2>/tmp/mtest.err & echo $!; wait $!; echo $? >> /tmp/{tmp['job_id']}/$!.exit_status)"
        print(cmd)
        popen_pipe = os.popen(cmd)
        print_file(cookies_file_path)
        pid = popen_pipe.readline().rstrip()
        print(f"PID: {pid}")
        tmp["pid"] = pid
        tmp["state_dict"] = dict(idx=i + 1, n_errors=n_errors, results_dir=results_dir, state="store")
        job_context = json.dumps(tmp)
        return None, job_context


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
    return job_context


@celery_app.task(name="success")
@celery_wf(celery_app, wf1, "success")
def wf1_complete_succesfully(job_context: str):
    """
    Just mark the Job as "completed succesfully"

    :param job_context:
    :return:
    """
    tmp = json.loads(job_context)
    tmp = change_status(tmp, "success")
    append_text("complete_successfully")
    job_context = json.dumps(tmp)
    append_text(f"cleanup:")
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
    tmp = change_status(tmp, "error")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    # los errores ya están en tmp o se descargaría en download
    #  para qué usaría la información en status?
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
    tmp = change_status(tmp, "cancel")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    job_executor.cancel_job(tmp['pid'])
    job_executor.remove_job_workspace(tmp)
    append_text(f"cancelled")
    return job_context
