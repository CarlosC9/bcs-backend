import unittest
from time import sleep

from biobarcoding.jobs import JobExecutorAtResourceFactory
import json
import requests
from biobarcoding.rest import bcs_api_base

"""
req = {
  "resource_id": "8fac3ce8-8796-445f-ac27-4baedadeff3b",
  "process_id": "c8df0c20-9cd5-499b-92d4-5fb35b5a369a",
  "process_params": {
    "parameters": {
      "ClustalW": {
        "darna": "PROTEIN"
      }
    },
    "data":  [
                    {
                      "step": "Input dataset",
                      "path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta",
                      "type": "fasta"
                    }
                  ]
  },
  "credentials": {
    "api_key": "fakekey"
  }
}



s = requests.Session()
s.put('http://localhost:5000/api/jobs/1?status=succes')
s.post('http://localhost:5000/api/jobs/', json = dict(params = req))
"""


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
    job_executor = JobExecutorAtResourceFactory().get(tmp["resource"]["jm_type"], tmp['resource'])
    job_executor.create_job_workspace(tmp['job_id'])
    job_context = json.dumps(tmp) # esto me puede dar como un output un id pero tmp me hace falta sabiendo el nombre
    # TODO comtemplar problema de conexión?
    # TODO Task "work"
    # TODO update Job status to "preparing_workspace"
    #  requests.put(job_context["endpoint_url"] + f"/api/jobs/{job_context['job_id']}/status", "preparing_workspace")
    return job_context


def wf1_transfer_data_to_resource(job_context: str):
    """
    Transfer data to the compute resource
    :param job_context:
    :return:
    """
    while(True):
        tmp = json.loads(job_context)
        filepackage_list = tmp["process"]['inputs']['data']
        tmp["process"]["workspace"] = "1cd8e2f6b131e891" #esto tiene que venir en el job_context

        job_executor = JobExecutorAtResourceFactory().get(tmp["resource"]["jm_type"], tmp['resource'])

        transfer_state = tmp.get("transfer_state")
        print(f"Transfer state: {transfer_state}")
        if transfer_state: # ya ha empezado la transferencia
            i = transfer_state["idx"]
            pid = transfer_state["pid"] # job_id
        else:
            i = 0
            pid = None

        # Ith transfer
        # transfer_at_i = tmp["process"]["inputs"]["parameters"]["upload_files"][i] if i < len(filepackage_list) else dict(
        #     local_path="", remote_path="")
        transfer_at_i = tmp["process"]['inputs']['data'][i] if i < len(filepackage_list) else dict(
            local_path="", remote_path="")


        local_path = transfer_at_i.get("path")
        remote_path = transfer_at_i.get("remote_path")  # Add workspace base?
        step = transfer_at_i.get("step")

        if i == len(filepackage_list):  # Transfer finished
            print("Transfer finished")
            del tmp["transfer_state"]
            job_context = json.dumps(tmp)
            return job_context
        elif job_executor.job_status(pid) == "running":  # Transfer is being executed
            # print("Transfer executing")
            sleep(5)
            print(job_context)
        elif job_executor.exists(local_path = local_path,
                                 remote_path = remote_path,
                                 step = step,
                                 workspace = str(tmp['job_id'])):  # File i has been transferred successfully
            print(f"File {i} transferred: : {local_path} -> {remote_path}. Moving to next")
            i += 1
            tmp["transfer_state"] = dict(idx=i, pid=None)
            job_context = json.dumps(tmp)
            print(job_context)
        else:  # Transfer file i
            print(f"Begin transfer {i}: {local_path} -> {remote_path}")
            pid = job_executor.upload_file(local_path=local_path,
                                            remote_path=remote_path,
                                            step=step,
                                            workspace=str(tmp['job_id']))
            # TODO yo puedo tener un error aquí
            tmp["transfer_state"] = dict(idx=i, pid=pid)
            print(tmp['transfer_state'])
            job_context = json.dumps(tmp)
            print(job_context)



def wf1_1(job_context):
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory()
    job_executor = job_executor.get(tmp["resource"]["jm_type"], tmp['resource'])
    inputs = tmp["process"]
    inv_id = job_executor.submit(str(tmp['job_id']), inputs)
    tmp['g_id'] = inv_id
    job_context = json.dumps(tmp)
    # append_text(outfile, f"submit. workspace: {tmp['g_id']}")
    return job_context


def wf1_2(job_context):
    tmp = json.loads(job_context)
    job = JobExecutorAtResourceFactory()
    job_executor = job.get(tmp["resource"]["jm_type"], tmp["resource"])
    status = job_executor.job_status(tmp["g_id"])
    if status == 'running' or 'ok':  # pasa siempre or running?
        return job_context
    if status == 'error':
        return 'error', job_context
    else:
        return 3, job_context


def wf1_3(job_context):
    tmp = json.loads(job_context)
    job = JobExecutorAtResourceFactory()
    job_executor = job.get(tmp["resource"]["jm_type"], tmp["resource"])
    status = job_executor.job_status(tmp["g_id"])
    if isinstance(status, dict):
        return 'error', job_context
    if status == 'ok':
        return job_context
    else:
        return 3, job_context


def wf1_errors(job_context):
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory()
    job_executor = job_executor.get(tmp["resource"]["jm_type"], tmp["resource"])
    status = job_executor.job_status(job_executor["g_id"])
    tmp['error'] = status
    job_context = json.dumps(tmp)
    return job_context


def wf_results(job_context):
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory()
    job_executor = job_executor.get(tmp["resource"]["jm_type"], tmp["resource"])
    r = job_executor.get_results(tmp["g_id"])
    tmp['results'] = r
    print(tmp['results'])
    job_context = json.dumps(tmp)
    return job_context


def wf_clean(job_context):
    tmp = json.loads(job_context)
    job = JobExecutorAtResourceFactory()
    job_executor = job.get(tmp["resource"]["jm_type"], tmp["resource"])
    job_executor.remove_job_workspace(str(tmp["job_id"]))
    del tmp['g_id']
    job_context = json.dumps(tmp)
    return job_context


def wf_cancel(job_context):
    tmp = json.loads(job_context)
    job_executor = JobExecutorAtResourceFactory()
    job_executor = job_executor.get(tmp["resource"]["jm_type"], tmp["resource"])
    job_executor.cancel_job(tmp['g_id'])
    return job_context


def workflow(job_context):
    # job_context = json.dumps(dict_input)
    return_0 = wf1_prepare_workspace(job_context)
    return_01 = wf1_transfer_data_to_resource(return_0)
    return_1 = wf1_1(return_01)
    if isinstance(return_1, tuple):
        print('error')
    else:
        return_2 = wf1_2(return_1)
        while isinstance(return_2, tuple):
            return_2 = wf1_2(return_1)
        return_3 = wf1_3(return_2)
        while isinstance(return_3, tuple):
            error, _ = return_3
            if error == 'error':
                return_3 = wf_cancel(return_2)
                return_4 = wf_clean(return_3)
                return return_4
            return_3 = wf1_3(return_2)
        if isinstance(return_3, str):
            return_4 = wf_results(return_3)
            return_5 = wf_clean(return_4)
            return return_5
        else:
            return_6 = wf1_errors(job_context)
            return return_6

class MyTestCase(unittest.TestCase):
    def test_celeryworkflow(self):

        #job_context = '{"endpoint_url": "", "process": {"inputs": {"parameters": {"ClustalW": {"darna": "PROTEIN"}}, "data": {"Input dataset": {"path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta", "type": "fasta"}}}, "name": "MSA ClustalW"}, "resource": {"name": "localhost - galaxy", "jm_type": "galaxy", "jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}}, "job_id": 52}'
        # job_context = '{"endpoint_url": "", "process": {"inputs": {"parameters": {"clustalw": {"dnarna": "DNA", "outform": "clustal", "out_order": "ALIGNED", "mode": "complete", "out_seqnos": "ON"}, "phyml": {"phylip_format": "", "nb_data_set": "1", "type_of_seq": "nt", "prop_invar": "e", "equi_freq": "m", "nbSubstCat": "4", "gamma": "e", "move": "NNI", "optimisationTopology": "tlr", "branchSupport": "-4", "numStartSeed": "0", "inputTree": "false", "tstv": "e", "model": "HKY85"}}, "data": {"Input dataset": {"path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta", "type": "fasta"}}}, "name": "ClustalW-PhyMl"}, "resource": {"name": "localhost - galaxy", "jm_type": "galaxy", "jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}}, "job_id": 8}'
        # NUEVO:
        job_context = '{"endpoint_url": "http//:localhost:5000/", "process": {"inputs": {"parameters": {"ClustalW": {"darna": "PROTEIN"}}, "data": [{"step": "Input dataset","path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta","type": "fasta"}]},"name": "MSA ClustalW"}, "status": "created","resource": {"name": "localhost - galaxy", "jm_type": "galaxy","jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}},"job_id": 60}'
        print(workflow(job_context))

if __name__ == '__main__':
    MyTestCase.test_celeryworkflow()
