import unittest
from time import sleep

from biobarcoding.jobs import JobExecutorAtResourceFactory
import json
import requests
from biobarcoding.rest import bcs_api_base
import requests
import os
import subprocess

from biobarcoding.rest.file_manager import FilesAPI
from biobarcoding.tasks.definitions import change_status
# req = {
#   "resource_id": "8fac3ce8-8796-445f-ac27-4baedadeff3b",
#   "process_id": "c8df0c20-9cd5-499b-92d4-5fb35b5a369a",
#   "process_params": {
#     "parameters": {
#       "ClustalW": {
#         "darna": "PROTEIN"
#       }
#     },
#     "data":  [
#                     {
#                       "step": "Input dataset",
#                       "path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta",
#                       "type": "fasta"
#                     }
#                   ]
#   },
#   "credentials": {
#     "api_key": "fakekey"
#   }
# }
#
#
#
# s = requests.Session()
# s.put('http://localhost:5000/api/jobs/1?status=succes')
# s.post('http://localhost:5000/api/jobs/', json = req)
# """
#
# """

TEST_DOWNLOAD_FOLDER = '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/inputs'
TEST_INPUTS_FOLDER = '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/inputs'
TEST_JOB_STATUS_DIR = '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/job_status'


def export(data,file) -> object:
    """
    Execute remote client script
    @param url: url of the get
    @param path: local path of the file gotten with curl
    @return: pid: PID of the executed script process
    """
    tmp_path = os.path.join(TEST_INPUTS_FOLDER, file)
    extension = data['type']
    bos_type = data['bo_type']
    selection_dict = {"filter": [{"feature_id": {"op": "in", "unary": data['selection']}}]}
    selection_json = json.dumps(selection_dict)
    endpoint = "http://localhost:5000"
    url = f"{endpoint}{bcs_api_base}/bos/{bos_type}.{extension}"
    curl = f"curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt -X GET -d \'{selection_json}\' -H \'Content-Type:application/json\' {url} -o {tmp_path}"
    cmd = f"(nohup bash -c &quot;{curl}&quot; >/tmp/mtest </dev/null 2>/tmp/mtest.err & echo $!; wait $!; echo $?  >> {TEST_JOB_STATUS_DIR}/$!.exit_status)"
    print(cmd)
    popen_pipe = os.popen(cmd)
    pid = popen_pipe.readline().rstrip()
    print(f"PID: {pid}")
    return pid

def export_status(pid):
    exit_status = "none"
    if pid:
        if os.path.isfile(f"{TEST_JOB_STATUS_DIR}/{pid}.exit_status"):
            with open(f"{TEST_JOB_STATUS_DIR}/{pid}.exit_status", "r") as f:
                exit_status = f.readline().strip()
                if exit_status.strip() == "0":
                    exit_status = "ok"
                else:
                    print(f"Error executing get with pid: {pid}. Exit status = {exit_status}")
                    exit_status = ""  # This means error
        else:
            exit_status = "running"

    return exit_status

def is_input_file(input_file):
    if os.path.exists(input_file): #and check that it is not an error file:
        return True
    else:
        return False



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
    # from saved JSON
    # curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt -X GET -d @"$TEST_FILES_PATH/seq_request.json" -H "Content-Type: application/json" http://localhost:5000/api/bos/sequences.fasta -o $TEST_FILES_PATH/test.fasta
    # from string
    # curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt -X GET -d '{"filter": [{"feature_id": {"op": "in", "unary": [1, 2]}}]}' -H "Content-Type:application/json"  http://localhost:5000/api/bos/sequences.fasta -o $TEST_FILES_PATH/test.fasta

    """
    Once the computation ends, transfer results from the resource to local
    :param job_context:
    :return:
    """
    MAX_ERRORS = 3
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    response = job_executor.check()
    print(response)
    transfer_state = tmp.get("transfer_state")
    # TODO acordar el path con Rafa
    base_path_to_tmp = TEST_INPUTS_FOLDER
    inputs_dir = os.path.join(base_path_to_tmp, str(tmp['job_id']))
    print(f"Transfer state: {transfer_state}")
    if transfer_state:  # ya ha empezado la transferencia
        i = transfer_state["idx"]
        n_errors = transfer_state["n_errors"]
    else:
        i = 0
        n_errors = 0
        # os.mkdir(inputs_dir)
        tmp["pid"] = None
        tmp["transfer_state"] = dict(idx=i, n_errors=n_errors, state="download", results_dir=inputs_dir)

    # Number of files to transfer
    n_files = len(tmp['process']['data'])

    if i < n_files and n_errors >= MAX_ERRORS:
        print(f"Transfer error: File {i + 1}")
        job_context = json.dumps(tmp)
        return 3, job_context
    if i == n_files:  # Transfer finished
        print("Transfer finished")
        del tmp["transfer_state"]
        job_context = json.dumps(tmp)
        return job_context
    elif export_status(tmp['pid']) == "running":  # Transfer is being executed
        print("Transfer executing")
        return 1, job_context
    elif export_status(tmp['pid']) == "":
        tmp["pid"] = None
        tmp["transfer_state"] = dict(idx=i, n_errors=n_errors + 1, state="download", results_dir=inputs_dir)
        job_context = json.dumps(tmp)
        return None, job_context
    elif is_input_file(os.path.join(TEST_INPUTS_FOLDER,job_executor.get_export_path(tmp))):
        print(f"File {i + 1} transferred -> Moving to next")
        i += 1
        tmp["pid"] = None
        tmp["transfer_state"] = dict(idx=i, n_errors=n_errors, state="download", results_dir=inputs_dir)
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i + 1}")
        file = job_executor.get_export_path(tmp)
        data = tmp['process']['data'][i]
        pid = export(data, file)
        tmp["pid"] = pid
        tmp["transfer_state"] = dict(idx=i, n_errors=0, state="download", results_dir=inputs_dir)
        job_context = json.dumps(tmp)
        return None, job_context



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
    job_context = json.dumps(tmp)
    return job_context


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
        pid = tmp["pid"]  # job_id
        n_errors = transfer_state["n_errors"]
    else:
        i = 0
        n_errors = 0
        tmp["pid"] = None
        tmp["transfer_state"] = dict(idx=i, n_errors=n_errors,state="upload")

    # Ith transfer
    files_list = job_executor.get_upload_files_list(tmp)
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
        del tmp['pid'] # TODO cambiar en definition
        del tmp["transfer_state"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.job_status(tmp) == "running":  # Transfer is being executed
        print("Transfer executing")
        sleep(5)
        return None, job_context
    elif job_executor.exists(tmp):  # File i has been transferred successfully
        print(f"File {i} transferred: {local_path} . Moving to next")
        i += 1
        tmp["pid"] = None
        tmp["transfer_state"] = dict(idx=i, n_errors=n_errors, state="upload")
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i}: {local_path}")
        tmp['ìd'] = job_executor.upload_file(tmp)
        # TODO yo puedo tener un error aquí
        tmp["transfer_state"] = dict(idx=i, n_errors=n_errors, state="upload")
        print(tmp['transfer_state'])
        job_context = json.dumps(tmp)
        return None, job_context


def wf1_submit(job_context: str):
    """
    Submit job to compute resource
    :param job_context:
    :return:
    """
    tmp = json.loads(job_context)
    tmp = change_status(tmp, "submit")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    pid = job_executor.submit(tmp) # TODO corregir en definition
    tmp['pid'] = pid
    job_context = json.dumps(tmp)
    print(f"submit. workspace: {tmp['pid']}")
    return job_context

def wf1_wait_until_execution_starts(job_context: str):
    """
     Wait for the job to start executing at the
     :param job_context:
     :return:
     """
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    status = job_executor.job_status(tmp)
    if isinstance(status, dict):  # todo esto sirve a los dos?
        print(f"wait_for_execution_end: status: {status}")
        tmp['process']['error'] = status
        job_context = json.dumps(tmp)
        return 'error', job_context
    elif status == 'running':
        print(f"wait_for_execution_end: status: {status}")
        return job_context
    else:
        print(f"wait_for_execution_end: status: {status}")
        return 3, job_context

def wf1_wait_for_execution_end(job_context: str):
    """
    Wait for the job to finish execution, knowing it is running
    When finished, it can end successfully or with an error.

    :param job_context:
    :return:
    """
    tmp = json.loads(job_context) # solo para test
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    status = job_executor.job_status(tmp)
    if isinstance(status, dict):
        print(f"wait_for_execution_end: status: {status}")
        return 'error', job_context
    elif status == 'ok':
        print(f"wait_for_execution_end: status: {status}")
        return job_context
    else:
        print(f"wait_for_execution_end: status: {status}")
        return 3, job_context


def wf1_transfer_data_from_resource(job_context: str):
    MAX_ERRORS = 3
    """
    Once the computation ends, transfer results from the resource to local
    :param job_context:
    :return:
    """
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    transfer_state = tmp.get("transfer_state")
    tmp = change_status(tmp,"transfer_data_from")
    # TODO acordar el path con Rafa
    base_path_to_results = '/tmp/'
    results_dir = os.path.join(base_path_to_results, str(tmp['job_id']))
    print(f"Transfer state: {transfer_state}")
    if transfer_state:  # ya ha empezado la transferencia
        i = transfer_state["idx"]
        n_errors = transfer_state["n_errors"]
    else:
        i = 0
        n_errors = 0
        # os.mkdir(results_dir)
        tmp["pid"] = None
        tmp["transfer_state"] = dict(idx=i, n_errors=n_errors, state="download", results_dir=results_dir)

    # Number of files to transfer
    n_files = len(job_executor.get_download_files_list(tmp))

    if i < n_files and n_errors >= MAX_ERRORS:
        print(f"Transfer error: File {i + 1}")
        job_context = json.dumps(tmp)
        return 3, job_context
    if i == n_files:  # Transfer finished
        print("Transfer finished")
        del tmp["transfer_state"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.job_status(tmp) == "running":  # Transfer is being executed
        print("Transfer executing")
        return 1, job_context
    elif job_executor.job_status(tmp) == "": #hay un error de descarga y lo vuelvo a intentar otra vez
        tmp["pid"] = None
        tmp["transfer_state"] = dict(idx=i, n_errors=n_errors + 1, state="download", results_dir=results_dir)
        job_context = json.dumps(tmp)
        return None, job_context
    elif job_executor.exists(tmp):  # File i has been transferred successfully
        print(f"File {i + 1} transferred -> Moving to next")
        # FilesAPI.put(local_path)
        i += 1
        tmp["pid"] = None
        tmp["transfer_state"] = dict(idx=i, n_errors=n_errors, state="download", results_dir=results_dir)
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i + 1}")
        pid = job_executor.download_file(tmp)
        tmp["pid"] = pid
        tmp["transfer_state"] = dict(idx=i, n_errors=0, state="download", results_dir=results_dir)
        job_context = json.dumps(tmp)
        return None, job_context


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
    print(f"cleanup:")
    return job_context


def wf1_complete_succesfully(job_context: str):
    """
      Just mark the Job as "completed succesfully"
      :param job_context:
      :return:
      """
    tmp = json.loads(job_context)
    tmp = change_status(tmp, "success")
    print("complete_successfully")
    job_context = json.dumps(tmp)
    print(f"cleanup:")
    return job_context


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
    print(f"error: {tmp['error']}")
    return job_context

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
    print(f"cancelled")
    return job_context

def workflow(job_context):
    # return_0 = wf1_prepare_workspace(job_context)
    # return_01 = (None, return_0)
    # while isinstance(return_01, tuple):
    #     _,job_context = return_01
    #     return_01 = wf1_transfer_data_to_resource(job_context)
    #     if isinstance(return_01,str):
    #         break
    #     if 'error' in return_01:
    #         raise Exception('error')
    # return_1 = wf1_submit(return_01)
    # return_2 = (None, return_1)
    # while isinstance(return_2, tuple):
    #     _, job_context = return_2
    #     return_2 = wf1_wait_until_execution_starts(job_context)  # problemas aquí
    #     if isinstance(return_2,str):
    #         break
    #     if 'error' in return_2:
    #         raise Exception('error')
    # return_3 = (None,return_2)
    # while isinstance(return_3, tuple):
    #     _,job_context = return_3
    #     return_3 = wf1_wait_for_execution_end(job_context)
    #     if isinstance(return_3,str):
    #         break
    #     if 'error' in return_3:
    #         raise Exception('error')
    return_4 = (None,job_context)
    while isinstance(return_4, tuple):
        _,job_context = return_4
        return_4 = wf1_transfer_data_from_resource(job_context)
        if isinstance(return_4,str):
            break
        if 'error' in return_4:
            raise Exception('error')
    return return_4



class MyTestCase(unittest.TestCase):
    job_context = {"endpoint_url": "http//:localhost:5000/",
                       "process":
                           {"inputs":
                                {"parameters": {"ClustalW": {"darna": "PROTEIN"}},
                                 # esta es la parte que sería la del export
                                 "data": [{"remote_name": "Input dataset",
                                           "file": "60/input_dataset.fasta",
                                           "type": "fasta"}]},
                            },
                       "name": "MSA ClustalW",
                       "results": [{"remote_name":"ClustalW on data 1: clustal", # esto sale de mi procesador clustal
                                    "file": "60/clustal.clustal",
                                    "type": "clustal"}, # ? esto es info tmb del proceso
                                   {"remote_name": "ClustalW on data 1: dnd",
                                    "file": "60/dnd.dnd",
                                    "type": "dnd"}
                                   # # ? esto es info tmb del proceso
                                   ],
                       "status": "created",
                       "resource": {"name": "localhost - galaxy", "jm_type": "galaxy","jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}},
                       "job_id": 60}

    job_context_id = {"endpoint_url": "http//:localhost:5000/",
                       "process":
                           {"inputs":
                                {"parameters": {"ClustalW": {"darna": "PROTEIN"}},
                                 # esta es la parte que sería la del export
                                 },
                            "data": [{"remote_name": "Input_dataset",
                                      "selection": [1, 2],
                                      "file": "60/input_dataset.fasta",
                                      "type": "fasta",
                                      "bo_type": "sequences"}]
                            },
                       "name": "MSA ClustalW",
                       "results": [{"remote_name":"ClustalW on data 1: clustal", # esto sale de mi procesador clustal
                                    "path": "out/",
                                    "file_ext": "clustal"}, # ? esto es info tmb del proceso
                                   {"remote_name": "ClustalW on data 1: dnd",
                                    "path": "out/dnd.dnd",
                                    "file_ext": "dnd"}
                                   # # ? esto es info tmb del proceso
                                   ],
                       "status": "created",
                       "resource": {"name": "localhost - galaxy", "jm_type": "galaxy","jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}},
                       "job_id": 60}

    job_context = json.dumps(job_context)
    job_context_id = json.dumps(job_context_id)
    def test_celeryworkflow(self):

        #job_context = '{"endpoint_url": "", "process": {"inputs": {"parameters": {"ClustalW": {"darna": "PROTEIN"}}, "data": {"Input dataset": {"path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta", "type": "fasta"}}}, "name": "MSA ClustalW"}, "resource": {"name": "localhost - galaxy", "jm_type": "galaxy", "jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}}, "job_id": 52}'
        # job_context = '{"endpoint_url": "", "process": {"inputs": {"parameters": {"clustalw": {"dnarna": "DNA", "outform": "clustal", "out_order": "ALIGNED", "mode": "complete", "out_seqnos": "ON"}, "phyml": {"phylip_format": "", "nb_data_set": "1", "type_of_seq": "nt", "prop_invar": "e", "equi_freq": "m", "nbSubstCat": "4", "gamma": "e", "move": "NNI", "optimisationTopology": "tlr", "branchSupport": "-4", "numStartSeed": "0", "inputTree": "false", "tstv": "e", "model": "HKY85"}}, "data": {"Input dataset": {"path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta", "type": "fasta"}}}, "name": "ClustalW-PhyMl"}, "resource": {"name": "localhost - galaxy", "jm_type": "galaxy", "jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}}, "job_id": 8}'
        # NUEVO:



        # job_context = '{"endpoint_url": "http//:localhost:5000/", ' \
        #               '"process": {"inputs": {"parameters": {"ClustalW": {"darna": "PROTEIN"}}, "data": [{"step": "Input dataset","path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta","type": "fasta"}]},"name": "MSA ClustalW"}, "status": "created",' \
        #               '"resource": ' \
        #               '{"name": "localhost - galaxy", "jm_type": "galaxy","jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}},' \
        #               '"job_id": 60, "output" : [{"local_path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/download", "name": }]}}'




        print(workflow(self.job_context))

    def test_export(self):
        return_4 = (None, self.job_context_id)
        while isinstance(return_4, tuple):
            _, job_context = return_4
            return_4 = wf1_export_to_supported_file_formats(job_context)
            if isinstance(return_4, str):
                break
            if 'error' in return_4:
                raise Exception('error')
        return return_4


if __name__ == '__main__':
    MyTestCase.test_celeryworkflow()
    MyTestCase.test_export()
