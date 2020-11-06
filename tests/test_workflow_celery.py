import unittest
from biobarcoding.jobs import JobExecutorAtResourceFactory
import json


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
    if status == 'running':  # pasa siempre or running?
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
    return_1 = wf1_1(job_context)
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

        job_context = '{"endpoint_url": "", "process": {"inputs": {"parameters": {"ClustalW": {"darna": "PROTEIN"}}, "data": {"Input dataset": {"path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta", "type": "fasta"}}}, "name": "MSA ClustalW"}, "resource": {"name": "localhost - galaxy", "jm_type": "galaxy", "jm_location": {"url": "http://localhost:8080/"}, "jm_credentials": {"api_key": "fakekey"}}, "job_id": 52}'


        # job_context_dict = {
        #     "job_id": str(52),
        #     "endpoint_url": "http://localhost:8080/",
        #     "resource": {"name": "beauvoir3",
        #                  "jm_type": "galaxy",
        #                  "jm_location": {"url": "http://localhost:8080/"},
        #                  "jm_credentials": {"api_key": "fakekey"}
        #                  },
        #     "process": {"name": "MSA ClustalW",
        #                 "inputs":
        #                     {"parameters":
        #                          {"ClustalW": {"darna": "PROTEIN"}
        #                           },
        #                      "data": {"Input dataset":
        #                                   {
        #                                       "path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta",
        #                                       "type": "fasta"
        #                                       }
        #                               }
        #                      }
        #                 }
        # }
        print(workflow(job_context))


if __name__ == '__main__':
    MyTestCase.test_celeryworkflow()
