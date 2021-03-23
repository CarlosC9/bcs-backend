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
        pid = transfer_state["pid"]  # job_id
        n_errors = transfer_state["n_errors"]
    else:
        i = 0
        pid = None
        n_errors = 0
        tmp["transfer_state"] = dict(idx=i, pid=None, n_errors=n_errors,state = "upload")

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
    elif job_executor.job_status(tmp) == "running":  # Transfer is being executed
        print("Transfer executing")
        sleep(5)
        return None, job_context
    elif job_executor.exists(tmp):  # File i has been transferred successfully
        print(f"File {i} transferred: {local_path} . Moving to next")
        i += 1
        tmp["transfer_state"] = dict(idx=i, pid=None, n_errors=n_errors, state = "upload")
        job_context = json.dumps(tmp)
        return None, job_context
    else:  # Transfer file i
        print(f"Begin transfer {i}: {local_path}")
        pid = job_executor.upload_file(tmp)
        # TODO yo puedo tener un error aquí
        tmp["transfer_state"] = dict(idx=i, pid=pid, n_errors=n_errors, state = "upload")
        print(tmp['transfer_state'])
        job_context = json.dumps(tmp)
        return None, job_context


def wf1_submit(job_context: str):
    """
    Submit job to compute resource
    :param job_context:
    :return:
    """
    tmp = json.loads(job_context) # solo para test
    # tmp = change_status(tmp, "submit")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    inputs = tmp["process"]
    pid = job_executor.submit(tmp)
    tmp['pid'] = pid
    job_context = json.dumps(tmp)
    return job_context

def wf1_wait_until_execution_starts(job_context: str):
    """
    Wait for the job to start executing at the
    :param job_context:
    :return:
    """
    tmp = json.loads(job_context) # solo para test
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    status = job_executor.job_status(tmp)
    if isinstance(status, dict):
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
    tmp = json.loads(job_context) # solo para test
    #TODO hacer esto general
    base_path_to_results = '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/download'
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    transfer_state = tmp.get("transfer_state")
    print(f"Transfer state: {transfer_state}")
    if transfer_state:  # ya ha empezado la transferencia
        i = transfer_state["idx"]
        pid = transfer_state["pid"]  # job_id
    else:
        i = 0
        n_errors = 0 # TODO esto al final?
        tmp["transfer_state"] = dict(idx=i, pid=None)
        # TODO va a existir el caso en el vaya a crearlo dos veces?
        try:
            os.mkdir(f"{base_path_to_results}/{tmp['job_id']}/")
        except FileExistsError:
            print("folder already existst")

    # Ith transfer
    # files_list = tmp["process"]["inputs"]["results"] # en mi caso esta es una lista de ids que tengo que pedir a mi workspace
    files_list = tmp["results"] # todo mejor ponerlo en tmp["process"]["results"] pero no con inputs
    transfer_at_i = files_list[i] if i < len(files_list) else dict(
        local_path="", remote_path="")

    '''The transfer_at_i keys are unique for each job type. The only field that is mandatory
    is the path that refers to the local path of the file to be transferred'''
    local_path = os.path.join(f"{base_path_to_results}/{tmp['job_id']}/", transfer_at_i.get("path"))

    if i == len(files_list):  # Transfer finished
        print("Transfer finished")
        del tmp["transfer_state"]
        job_context = json.dumps(tmp)
        return job_context
    elif job_executor.job_status(tmp) == "running":  # Transfer is being executed
        print("Transfer executing")
        sleep(5)
        print(job_context)
    elif os.path.isfile(local_path):  # File i has been transferred successfully
        print(f"File {i} transferred: {local_path} -> Moving to next")
        FilesAPI.put(local_path)
        i += 1
        tmp["transfer_state"] = dict(idx=i, pid=None)
        job_context = json.dumps(tmp)
        print(job_context)
    else:  # Transfer file i
        print(f"Begin transfer {i}: {local_path}")
        pid = job_executor.download_file(tmp)
        # TODO yo puedo tener un error aquí
        tmp["transfer_state"] = dict(idx=i, pid=pid)
        print(tmp['transfer_state'])
        job_context = json.dumps(tmp)
        return None, job_context


def wf1_cleanup_workspace(job_context: str):
    """
    Delete workspace at the remote resource

    :param job_context:
    :return:
    """

    tmp = json.loads(json.loads(job_context)) # solo para test
    tmp = change_status(tmp, "cleanup")
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    job_executor.remove_job_workspace(str(tmp["job_id"]))
    del tmp['pid']
    job_context = json.dumps(tmp)
    # append_text(f"cleanup:")


def wf1_complete_succesfully(job_context: str):
    """
    Just mark the Job as "completed succesfully"

    :param job_context:
    :return:
    """
    # TODO CHANGE STATE
    # append_text("complete_successfully")
    return job_context


def wf1_completed_error(job_context: str):
    """
    Mark the Job as "completed with error"

    :param job_context:
    :return:
    """

    tmp = json.loads(json.loads(job_context)) # solo para test
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    status = job_executor.job_status(job_executor["pid"])
    tmp['error'] = status
    job_executor.remove_job_workspace(str(tmp['job_id']))
    job_context = json.dumps(tmp)
    # append_text(f"error: {tmp['error']}")
    return job_context

def wf1_cancelled(job_context: str):
    """
    Mark the Job as "cancelled"

    :param job_context:
    :return:
    """
    tmp = json.loads(json.loads(job_context)) # solo para test
    job_executor = JobExecutorAtResourceFactory().get(tmp)
    job_executor.cancel_job(tmp['pid'])
    job_executor.remove_job_workspace(tmp)
    # append_text(f"cancelled")
    return job_context

def workflow(job_context):
    return_0 = wf1_prepare_workspace(job_context)
    return_01 = (None, return_0)
    while isinstance(return_01, tuple):
        _,job_context = return_01
        return_01 = wf1_transfer_data_to_resource(job_context)
        if isinstance(return_01,str):
            break
        if 'error' in return_01:
            raise Exception('error')
    return_1 = wf1_submit(return_01)
    return_2 = (None, return_1)
    while isinstance(return_2, tuple):
        _, job_context = return_2
        return_2 = wf1_wait_until_execution_starts(job_context)  # problemas aquí
        if isinstance(return_2,str):
            break
        if 'error' in return_2:
            raise Exception('error')
    return_3 = (None,return_2)
    while isinstance(return_3, tuple):
        _,job_context = return_3
        return_3 = wf1_wait_for_execution_end(job_context)
        if isinstance(return_3,str):
            break
        if 'error' in return_3:
            raise Exception('error')
    return_4 = (None,return_3)
    while isinstance(return_4, tuple):
        _,job_context = return_4
        return_4 = wf1_transfer_data_from_resource(return_3)
        if isinstance(return_4,str):
            break
        if 'error' in return_4:
            raise Exception('error')
    return return_4



class MyTestCase(unittest.TestCase):
    def test_celeryworkflow(self):

        #job_context = '{"endpoint_url": "", "process": {"inputs": {"parameters": {"ClustalW": {"darna": "PROTEIN"}}, "data": {"Input dataset": {"path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta", "type": "fasta"}}}, "name": "MSA ClustalW"}, "resource": {"name": "localhost - galaxy", "jm_type": "galaxy", "jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}}, "job_id": 52}'
        # job_context = '{"endpoint_url": "", "process": {"inputs": {"parameters": {"clustalw": {"dnarna": "DNA", "outform": "clustal", "out_order": "ALIGNED", "mode": "complete", "out_seqnos": "ON"}, "phyml": {"phylip_format": "", "nb_data_set": "1", "type_of_seq": "nt", "prop_invar": "e", "equi_freq": "m", "nbSubstCat": "4", "gamma": "e", "move": "NNI", "optimisationTopology": "tlr", "branchSupport": "-4", "numStartSeed": "0", "inputTree": "false", "tstv": "e", "model": "HKY85"}}, "data": {"Input dataset": {"path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta", "type": "fasta"}}}, "name": "ClustalW-PhyMl"}, "resource": {"name": "localhost - galaxy", "jm_type": "galaxy", "jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}}, "job_id": 8}'
        # NUEVO:



        job_context = '{"endpoint_url": "http//:localhost:5000/", ' \
                      '"process": {"inputs": {"parameters": {"ClustalW": {"darna": "PROTEIN"}}, "data": [{"step": "Input dataset","path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta","type": "fasta"}]},"name": "MSA ClustalW"}, "status": "created",' \
                      '"resource": ' \
                      '{"name": "localhost - galaxy", "jm_type": "galaxy","jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}},' \
                      '"job_id": 60, "output" : [{"local_path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/download", "name": }]}}'

        job_context = {"endpoint_url": "http//:localhost:5000/",
                       "process":
                           {"inputs":
                                {"parameters": {"ClustalW": {"darna": "PROTEIN"}},
                                 # esta es la parte que sería la del export
                                 "data": [{"remote_name": "Input dataset",
                                           "path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta",
                                           "type": "fasta"}]},
                            },
                       "name": "MSA ClustalW",
                       "results": [{"remote_path":"clustal", # esto sale de mi procesador clustal
                                    "path": "/out/clustal.clustal",
                                    "file_ext": "clustal"}, # ? esto es info tmb del proceso
                                   {"remote_name": "dnd",
                                    "path": "/out/dnd.dnd",
                                    "file_ext": "dnd"} # # ? esto es info tmb del proceso
                                   ],
                       "status": "created",
                       "resource": {"name": "localhost - galaxy", "jm_type": "galaxy","jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}},
                       "job_id": 60}
        job_context = json.dumps(job_context)

        print(workflow(job_context))

if __name__ == '__main__':
    MyTestCase.test_celeryworkflow()
