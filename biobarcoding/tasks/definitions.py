import json
import os
import re
import subprocess
from biobarcoding.rest import bcs_api_base
from biobarcoding.tasks import celery_app
from biobarcoding.common.decorators import celery_wf
from biobarcoding.jobs import JobExecutorAtResourceFactory
from biobarcoding.common import ROOT

"""
Celery tasks CANNOT be debugged in Celery!! (the code is run in a separate process; 
of course they can be debugged "off-line")
"""

MAX_ATTEMPTS = 3
CELERY_LOG = ROOT + "/tests/data_test/celery_log.txt"


# Send messages
# Refresh queues of jobs (at different computing resources)

def change_status(tmp, status: str):
    endpoint_url = os.getenv("ENDPOINT_URL")
    cookies_file_path = os.getenv("COOKIES_FILE_PATH")
    job_id = tmp["job_id"]
    url = f"{endpoint_url}{bcs_api_base}/jobs/{job_id}"
    if tmp["status"] != status:
        api_login()
        status_request = json.dumps(dict(status=status))
        cmd = ["curl", "--cookie-jar", cookies_file_path, "--cookie", cookies_file_path, "-H",
               "Content-Type: application/json", "-XPUT", "--data-binary", status_request, url]
        print(subprocess.run(cmd))
        tmp["status"] = status
        return tmp
    else:
        return tmp


def write_to_file(filename, s):
    with open(filename, "a+") as f:
        f.write(s + "\r\n")


def create_file(filename):
    if not os.path.exists(filename):
        open(filename, "x")


def api_login():
    endpoint_url = os.getenv("ENDPOINT_URL")
    cookies_file_path = os.getenv("COOKIES_FILE_PATH")
    url = f"{endpoint_url}{bcs_api_base}/authn?user=test_user"
    cmd = ["curl", "--cookie-jar", cookies_file_path, "-X", "PUT", url]
    subprocess.run(cmd)


def check_file_is_stored_in_backend(filename, job_id):
    endpoint_url = os.getenv("ENDPOINT_URL")
    cookies_file_path = os.getenv("COOKIES_FILE_PATH")
    api_login()
    cmd = ["curl", "--cookie-jar", cookies_file_path, "--cookie",
           cookies_file_path, f"{endpoint_url}{bcs_api_base}/files/jobs/{job_id}/{filename}"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    process_return_dict = json.loads(proc.stdout)
    if process_return_dict.get("content"):
        return process_return_dict["content"].get("size") != 0
    else:
        return False


def put_logs_in_job_context_for_store(job_context, local_workspace):
    for f in os.listdir(local_workspace):
        if re.search(r"bcs\..+\.log", f):
            job_context["process"]["inputs"]["parameters"]["result_files"].append(
                {
                    "remote_name": "",
                    "file": f,
                    "type": "stdout" if "stdout" in f else "stderr"
                }
            )
    return job_context


def clean_failed_results(result_files, local_workspace):
    clean_results = []
    for f in result_files:
        if os.path.exists(os.path.join(local_workspace, f["file"])):
            clean_results.append(f)
    return clean_results


def clean_exit_status_from_local_workspase(local_workspace):
    for f in os.listdir(local_workspace):
        if re.search(r"\d+\.exit_status", f):
            os.remove(os.path.join(local_workspace, f))


def check_exit_status_from_local_workspace_is_cleaned(local_workspace):
    cleaned = True
    for f in os.listdir(local_workspace):
        if re.search(r"\d+\.exit_status", f):
            cleaned = False
    return cleaned


# TODO: Hay que prepararlo para las colecciones y los ficheros que suba el usuario de manera
# que se puedan concatenar si se refieren al mismo fichero.
# Mirar: https://stackoverflow.com/questions/40359012/how-to-append-a-file-with-the-existing-one-using-curl

def export(file_dict, job_executor) -> object:
    """
    Execute remote client script
    @param file_dict: Dictionary of the file to be exported
    @param job_executor: Needed to know local_workspace and log names
    @return: pid: PID of the executed script process
    """
    tmp_path = os.path.join(job_executor.local_workspace, file_dict["file"])
    extension = file_dict['type']
    bos_type = file_dict['bo_type']
    selection_dict = {"filter": [{"feature_id": {"op": "in", "unary": file_dict['selection']}}]}
    selection_json = json.dumps(selection_dict)
    endpoint = os.getenv("ENDPOINT_URL")
    cookies_file_path = os.getenv("COOKIES_FILE_PATH")
    api_login()
    url = f"{endpoint}{bcs_api_base}/bos/{bos_type}.{extension}"
    curl = f"curl --cookie-jar {cookies_file_path} --cookie {cookies_file_path} -X GET -d \'{selection_json}\' -H \'Content-Type:application/json\' {url} -o {tmp_path}"
    cmd = f"(nohup bash -c &quot;{curl}&quot; >{job_executor.log_filenames_dict['export_stdout']} </dev/null 2>{job_executor.log_filenames_dict['export_stderr']} & echo $!; wait $!; echo $?  >> {job_executor.local_workspace}/$!.exit_status)"
    print(cmd)
    popen_pipe = os.popen(cmd)
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


# To test:
# - Run "bcs-backend" in the IDE
# - Open three shells in the "bcs-backend" directory:
#   - ./start_local_dev_services.sh
#   - python3 biobarcoding/rest/main.py
#   - celery -A biobarcoding.tasks.definitions flower [OPTIONAL: to see tasks, http://localhost:5555]
#   - curl -i -XPOST http://localhost:5000/api/jobs/ --data-urlencode "{}" --> (should return immediately)
#     - (to see messages generated by tasks of this workflow) tail -f /home/rnebot/Downloads/borrame/log.txt
#   - curl -i -XPOST http://localhost:5000/api/jobs/ -H "Content-Type: application/json" -d @"/home/daniel/Documentos/GIT/bcs-backend/tests/request_transfer.json"
#   - curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt -XPOST http://localhost:5000/api/jobs/ -H "Content-Type: application/json" -d @"/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/new_galaxy_request.json"
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
    tmp = json.loads(job_context)
    tmp = change_status(tmp, "preparing_workspace")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    state_dict = tmp.get("state_dict")

    print(f"Prepare state: {state_dict}")

    if state_dict:  # ya ha empezado el prepare
        n_attempts = state_dict["n_attempts"]
    else:
        n_attempts = 0
        create_file(job_executor.log_filenames_dict["prepare_stdout"])
        create_file(job_executor.log_filenames_dict["prepare_stderr"])
        tmp["state_dict"] = dict(n_attempts=n_attempts, state="prepare")

    if not job_executor.check():
        error_str = "Connection to the server has been lost"
        tmp["error"] = error_str
        write_to_file(job_executor.log_filenames_dict["prepare_stderr"], error_str)
        job_context = json.dumps(tmp)
        return "error", job_context
    elif n_attempts >= MAX_ATTEMPTS:
        error_str = f"It was impossible to prepare the working directory of Job {tmp['job_id']}." + \
                    f"Maybe it could be due to some disk space or credentials issue."
        tmp["error"] = error_str
        write_to_file(job_executor.log_filenames_dict["prepare_stderr"], error_str)
        job_context = json.dumps(tmp)
        return "error", job_context
    elif job_executor.job_workspace_exists(str(tmp["job_id"])):
        write_to_file(job_executor.log_filenames_dict["prepare_stdout"],
                      f"Job {tmp['job_id']} workspace prepared")
        del tmp["state_dict"]
        job_context = json.dumps(tmp)
        return job_context
    else:
        write_to_file(job_executor.log_filenames_dict["prepare_stdout"],
                      f"Preparing Job {tmp['job_id']} workspace: Attempt: {n_attempts + 1}")
        tmp["state_dict"] = dict(n_attempts=n_attempts + 1, state="prepare")
        job_executor.create_job_workspace(str(tmp["job_id"]))
        job_context = json.dumps(tmp)
        return None, job_context


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
    tmp = change_status(tmp, "export")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    state_dict = tmp.get("state_dict")

    print(f"Export state: {state_dict}")

    if state_dict:  # ya ha empezado la transferencia
        i = state_dict["idx"]
        n_attempts = state_dict["n_attempts"]
    else:
        create_file(job_executor.log_filenames_dict["export_stdout"])
        create_file(job_executor.log_filenames_dict["export_stderr"])
        i = 0
        n_attempts = 0
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_attempts=n_attempts, state="export")

    # Files to transfer
    files = tmp['process']['inputs']['data']
    file_dict = files[i] if i < len(files) else {"file": "", "remote_name": "", "type": ""}

    if i < len(files) and n_attempts >= MAX_ATTEMPTS:
        error_str = f"Export error: File {i + 1} {file_dict['file']}"
        write_to_file(job_executor.log_filenames_dict["export_stderr"], error_str)
        tmp['error'] = error_str
        job_context = json.dumps(tmp)
        return "error", job_context
    if i == len(files):  # Transfer finished
        print("Export finished")
        write_to_file(job_executor.log_filenames_dict["export_stdout"], "Export finished successfully")
        del tmp["state_dict"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.step_status(tmp) == "running":  # Transfer is being executed
        print("Export executing")
        return 1, job_context
    elif job_executor.step_status(tmp) == "":
        tmp["pid"] = None
        job_context = json.dumps(tmp)
        return None, job_context
    elif is_bos_file(os.path.join(job_executor.local_workspace, file_dict["file"])):
        print(f"File {file_dict['file']} transferred -> Moving to next")
        write_to_file(job_executor.log_filenames_dict["export_stdout"],
                      f"File {file_dict['file']} exported -> Moving to next")
        tmp['process']['inputs']['data'][i]['file'] = os.path.join(job_executor.local_workspace, file_dict["file"])
        del tmp['process']['inputs']['data'][i]['selection']
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i + 1, n_attempts=0, state="export")
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin export {i + 1}: {file_dict['file']}. Attempt: {n_attempts + 1}")
        write_to_file(job_executor.log_filenames_dict["export_stdout"],
                      f"Begin export {i + 1}: {file_dict['file']}. Attempt: {n_attempts + 1}")
        pid = export(file_dict, job_executor)
        tmp["pid"] = pid
        tmp["state_dict"] = dict(idx=i, n_attempts=n_attempts + 1, state="export")
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
    state_dict = tmp.get("state_dict")
    print(f"Transfer state: {state_dict}")
    if state_dict:  # ya ha empezado la transferencia
        i = state_dict["idx"]
        n_attempts = state_dict["n_attempts"]
    else:
        create_file(job_executor.log_filenames_dict["upload_stdout"])
        create_file(job_executor.log_filenames_dict["upload_stderr"])
        i = 0
        n_attempts = 0
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_attempts=n_attempts, state="upload")

    # Files to transfer
    files = job_executor.get_upload_files_list(tmp)
    file_dict = files[i] if i < len(files) else {"file": "", "remote_name": "", "type": ""}

    if not job_executor.check():
        error_str = "Connection to the server has been lost"
        print(error_str)
        write_to_file(job_executor.log_filenames_dict["upload_stderr"], error_str)
        tmp['error'] = error_str
        job_context = json.dumps(tmp)
        return "error", job_context
    elif i < len(files) and n_attempts >= MAX_ATTEMPTS:
        error_str = f"Transfer of file {os.path.basename(file_dict['file'])} to resource failed."
        write_to_file(job_executor.log_filenames_dict["upload_stderr"], error_str)
        print(error_str)
        tmp['error'] = error_str
        job_context = json.dumps(tmp)
        return "error", job_context
    elif i == len(files):  # Transfer finished
        print("Transfer finished")
        write_to_file(job_executor.log_filenames_dict["upload_stdout"],
                      "Upload finished successfully")
        del tmp["state_dict"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.step_status(tmp) == "running":  # Transfer is being executed
        print("Transfer executing")
        return 1, job_context
    elif job_executor.step_status(tmp) == "":  # Transfer error
        print(f"Transfer error: File {i + 1}")
        tmp["pid"] = None
        job_context = json.dumps(tmp)
        return None, job_context
    elif job_executor.exists(tmp):  # File i has been transferred successfully
        print(f"File {file_dict['file']} transferred -> Moving to next. Attempt: {n_attempts}")
        write_to_file(job_executor.log_filenames_dict["upload_stdout"],
                      f"File {os.path.basename(file_dict['file'])} transferred -> Moving to next. Attempt: {n_attempts + 1}")
        i += 1
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_attempts=0, state="upload")
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i + 1}: {file_dict['file']}")
        write_to_file(job_executor.log_filenames_dict["upload_stdout"],
                      f"Begin transfer {i + 1}: {file_dict['file']}")
        pid = job_executor.upload_file(tmp)
        tmp["pid"] = pid
        tmp["state_dict"] = dict(idx=i, n_attempts=n_attempts + 1, state="upload")
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
    create_file(job_executor.log_filenames_dict["submit_stdout"])
    create_file(job_executor.log_filenames_dict["submit_stderr"])
    pid = job_executor.submit(tmp["process"])
    tmp['pid'] = pid
    tmp['state_dict'] = dict(state="submit", substep="submit")
    job_context = json.dumps(tmp)
    write_to_file(CELERY_LOG, f"Submit job with PID: {tmp['pid']}")
    write_to_file(job_executor.log_filenames_dict["submit_stdout"],
                  f"Submit job with PID: {tmp['pid']}")
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
    tmp = change_status(tmp, "wait_until_execution_starts")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    tmp['state_dict'] = dict(state="submit", substep="wait_until_execution_starts")
    if not job_executor.check():
        error_str = "Connection to the server has been lost"
        write_to_file(job_executor.log_filenames_dict["submit_stderr"], error_str)
        print(error_str)
        tmp['error'] = error_str
        job_context = json.dumps(tmp)
        return "error", job_context
    status = job_executor.step_status(tmp)
    write_to_file(CELERY_LOG, f"wait_until_execution_starts: status: {status}")
    # TODO: el step_status debería de sacar outputs uniformes. Se podría escribir ese
    # status en el log en el job_executor?
    if isinstance(status, dict) or status == "":
        if status != "":
            tmp['process']['error'] = status
            error_str = f"The process failed to start with status: {status}."
        else:
            error_str = "The process failed to start."
        print(error_str)
        write_to_file(job_executor.log_filenames_dict["submit_stderr"], error_str)
        tmp['error'] = error_str
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
    tmp = change_status(tmp, "wait_for_execution_end")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    tmp['state_dict'] = dict(state="submit", substep="wait_for_execution_end")
    status = job_executor.step_status(tmp)
    write_to_file(CELERY_LOG, f"wait_for_execution_end: status: {status}")
    if not job_executor.check():
        error_str = "Connection to the server has been lost"
        print(error_str)
        write_to_file(job_executor.log_filenames_dict["submit_stderr"], error_str)
        tmp['error'] = error_str
        job_context = json.dumps(tmp)
        return "error", job_context
    elif isinstance(status, dict):
        tmp['process']['error'] = status
        write_to_file(job_executor.log_filenames_dict["submit_stderr"], status)
        job_context = json.dumps(tmp)
        return 'error', job_context
    elif status == 'ok':
        del tmp["state_dict"]
        write_to_file(job_executor.log_filenames_dict["submit_stdout"],
                      "Job execution finished successfully.")
        job_context = json.dumps(tmp)
        return job_context
    elif status == "":
        job_executor.write_submit_logs()
        write_to_file(job_executor.log_filenames_dict["submit_stdout"],
                      "Job execution finished in error.")
        return 'error', job_context
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
    tmp = change_status(tmp, "transfer_data_from_resource")
    job_executor = JobExecutorAtResourceFactory().get(tmp)

    state_dict = tmp.get("state_dict")
    print(f"Transfer state: {state_dict}")
    if state_dict:  # ya ha empezado la transferencia
        i = state_dict["idx"]
        n_attempts = state_dict["n_attempts"]
    else:
        create_file(job_executor.log_filenames_dict["download_stdout"])
        create_file(job_executor.log_filenames_dict["download_stderr"])
        i = 0
        n_attempts = 0
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_attempts=n_attempts, state="download")

    # Files to transfer
    files = job_executor.get_download_files_list(tmp)
    file_dict = files[i] if i < len(files) else {"file": "", "remote_name": "", "type": ""}

    if not job_executor.check():
        error_str = "Connection to the server has been lost"
        print(error_str)
        write_to_file(job_executor.log_filenames_dict["download_stderr"], error_str)
        tmp['error'] = error_str
        job_context = json.dumps(tmp)
        return "error", job_context
    elif i < len(files) and n_attempts >= MAX_ATTEMPTS:
        error_str = f"Transfer from resource error: File {file_dict['file']}"
        write_to_file(job_executor.log_filenames_dict["download_stderr"], error_str)
        print(error_str)
        tmp['error'] = error_str
        job_context = json.dumps(tmp)
        return "error", job_context
    if i == len(files):  # Transfer finished
        print("Transfer finished")
        write_to_file(job_executor.log_filenames_dict["download_stdout"], "Download finished successfully")
        del tmp["state_dict"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.step_status(tmp) == "running":  # Transfer is being executed
        print("Transfer executing")
        return 1, job_context
    elif job_executor.step_status(tmp) == "":
        tmp["pid"] = None
        job_context = json.dumps(tmp)
        return None, job_context
    elif job_executor.exists(tmp):  # File i has been transferred successfully
        print(f"File {file_dict['file']} transferred -> Moving to next")
        write_to_file(job_executor.log_filenames_dict["download_stdout"],
                      f"File {file_dict['file']} transferred -> Moving to next. Attempts: {n_attempts}")
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i + 1, n_attempts=0, state="download")
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i + 1}: {file_dict['file']}. Attempt: {n_attempts + 1}")
        write_to_file(job_executor.log_filenames_dict["download_stdout"],
                      f"Begin transfer {i + 1}: {file_dict['file']}. Attempt: {n_attempts + 1}")
        pid = job_executor.download_file(tmp)
        tmp["pid"] = pid
        tmp["state_dict"] = dict(idx=i, n_attempts=n_attempts + 1, state="download")
        job_context = json.dumps(tmp)
        return None, job_context


@celery_app.task(name="store_result_in_backend")
@celery_wf(celery_app, wf1, "store_result_in_backend")
def wf1_store_result_in_backend(job_context: str):
    tmp = json.loads(job_context)
    tmp = change_status(tmp, "cleanup")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    state_dict = tmp.get("state_dict")

    print(f"Store result in backend state: {state_dict}")
    if state_dict:  # ya ha empezado la transferencia
        i = state_dict["idx"]
        n_attempts = state_dict["n_attempts"]
    else:
        tmp = put_logs_in_job_context_for_store(tmp, job_executor.local_workspace)
        create_file(job_executor.log_filenames_dict["store_stdout"])
        create_file(job_executor.log_filenames_dict["store_stderr"])
        i = 0
        n_attempts = 0
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i, n_attempts=n_attempts, state="store")

    # Files to transfer
    files = job_executor.get_download_files_list(tmp)
    file_dict = files[i] if i < len(files) else {"file": "", "remote_name": "", "type": ""}

    if i < len(files) and n_attempts >= MAX_ATTEMPTS:
        error_str = f"Store result in backend error: File {file_dict['file']}"
        write_to_file(job_executor.log_filenames_dict["store_stderr"], error_str)
        print(error_str)
        tmp['error'] = error_str
        job_context = json.dumps(tmp)
        return "error", job_context
    elif i == len(files):  # Transfer finished
        if files[-1]["file"] != "bcs.store.stderr.log":
            tmp["process"]["inputs"]["parameters"]["result_files"] += [
                {
                    "remote_name": "",
                    "file": "bcs.store.stdout.log",
                    "type": "stdout"
                },
                {
                    "remote_name": "",
                    "file": "bcs.store.stderr.log",
                    "type": "stderr"
                },
            ]
            job_context = json.dumps(tmp)
            return None, job_context

        print("Store result in backend finished")
        write_to_file(job_executor.log_filenames_dict["store_stdout"], "Store result in backend finished successfully")
        del tmp["state_dict"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.step_status(tmp) == "running":  # Transfer is being executed
        print("Store result in backend executing")
        return 1, job_context
    elif job_executor.step_status(tmp) == "":
        tmp["pid"] = None
        job_context = json.dumps(tmp)
        return None, job_context
    elif check_file_is_stored_in_backend(file_dict["file"], tmp["job_id"]):  # File i has been transferred successfully
        print(f"File {file_dict['file']} stored -> Moving to next")
        write_to_file(job_executor.log_filenames_dict["store_stdout"],
                      f"File {file_dict['file']} stored -> Moving to next. Attempts: {n_attempts}")
        tmp["pid"] = None
        tmp["state_dict"] = dict(idx=i + 1, n_attempts=0, state="store")
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Beginning to store result in backend of file {file_dict['file']}. Attempt: {n_attempts + 1}")
        write_to_file(job_executor.log_filenames_dict["store_stdout"],
                      f"Beginning to store result in backend of file {file_dict['file']}. Attempt: {n_attempts + 1}")
        cookies_file_path = os.getenv("COOKIES_FILE_PATH")
        endpoint_url = os.getenv("ENDPOINT_URL")
        local_path = os.path.join(job_executor.local_workspace, file_dict["file"])
        api_login()
        curl_cmd = f"curl -s --cookie-jar {cookies_file_path} --cookie {cookies_file_path} -H \"Content-Type: application/x-fasta\" -XPUT --data-binary @\"{local_path}\" \"{endpoint_url}{bcs_api_base}/files/jobs/{str(tmp['job_id'])}/{file_dict['file']}.content\""
        cmd = f"(nohup bash -c \'{curl_cmd} \' >/tmp/mtest </dev/null 2>/tmp/mtest.err & echo $!; wait $!; echo $? >> /tmp/{tmp['job_id']}/$!.exit_status)"
        print(cmd)
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().rstrip()
        print(f"PID: {pid}")
        tmp["pid"] = pid
        tmp["state_dict"] = dict(idx=i, n_attempts=n_attempts + 1, state="store")
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
    state_dict = tmp.get("state_dict")

    print(f"Cleanup state: {state_dict}")

    if state_dict:  # ya ha empezado el prepare
        n_attempts = state_dict["n_attempts"]
    else:
        n_attempts = 0
        create_file(job_executor.log_filenames_dict["cleanup_stdout"])
        create_file(job_executor.log_filenames_dict["cleanup_stderr"])
        tmp["state_dict"] = dict(n_attempts=n_attempts)

    if not job_executor.check():
        error_str = "Connection to the server has been lost"
        write_to_file(job_executor.log_filenames_dict["cleanup_stderr"], error_str)
        print(error_str)
        tmp['error'] = error_str
        job_context = json.dumps(tmp)
        return "error", job_context
    elif n_attempts >= MAX_ATTEMPTS:
        error_str = f"It was impossible to prepare the working directory of Job {tmp['job_id']}." + \
                    f"Maybe it could be due to some disk space or credentials issue."
        tmp["error"] = error_str
        write_to_file(job_executor.log_filenames_dict["cleanup_stderr"], error_str)
        job_context = json.dumps(tmp)
        return "error", job_context
    elif (not job_executor.job_workspace_exists(str(tmp["job_id"])) and
          check_exit_status_from_local_workspace_is_cleaned(job_executor.local_workspace)):
        write_to_file(job_executor.log_filenames_dict["cleanup_stdout"],
                      f"Job {tmp['job_id']} workspace cleaned up")
        del tmp["state_dict"]
        write_to_file(CELERY_LOG, f"cleanup:")
        job_context = json.dumps(tmp)
        return job_context
    else:
        write_to_file(job_executor.log_filenames_dict["cleanup_stdout"],
                      f"Cleaning Job {tmp['job_id']} workspace: Attempt: {n_attempts + 1}")
        job_executor.remove_job_workspace(str(tmp["job_id"]))
        clean_exit_status_from_local_workspase(job_executor.local_workspace)
        tmp["state_dict"] = dict(n_attempts=n_attempts + 1)
        job_context = json.dumps(tmp)
        return None, job_context


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
    write_to_file(CELERY_LOG, "complete_successfully")
    job_context = json.dumps(tmp)
    write_to_file(CELERY_LOG, f"success:")
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
    write_to_file(CELERY_LOG, f"error: {tmp['error']}")
    previous_state = tmp["state_dict"]["state"]
    if previous_state == "store" or previous_state == "cleanup":
        return job_context

    job_executor = JobExecutorAtResourceFactory().get(tmp)
    result_files = tmp["process"]["inputs"]["parameters"]["result_files"]
    result_files = clean_failed_results(result_files, job_executor.local_workspace)
    tmp["process"]["inputs"]["parameters"]["result_files"] = result_files
    del tmp["state_dict"]#needed to store to work
    job_context = json.dumps(tmp)

    return "wf1_store_result_in_backend", job_context


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
    write_to_file(CELERY_LOG, f"cancelled")
    return job_context
