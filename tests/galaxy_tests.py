import unittest
from pathlib import Path
import os

from bioblend import galaxy
from biobarcoding.jobs.galaxy_resource import *


class MyTestCase(unittest.TestCase):
    def test_upload_file(self):
        insfile = 'data_test/parsec_creds.yaml'
        insname = 'beauvoir3'
        ins = galaxy_instance(insfile, name = insname)
        gi = login(ins['key'],url=ins['url'])
        history = gi.histories.create_history(name="test_upload_file history")
        fn =  Path('data_test/ls_orchid.fasta')
        file_name = "test1"
        d = gi.tools.upload_file(
            fn,
            history_id=history["id"],
            file_name=file_name,
            dbkey="?",
            file_type="fasta")
        self.assertNotEqual(d, None, "should be something here")

    def test_inputs_files(self):
        insfile = 'data_test/parsec_creds.yaml'
        insname = 'local'
        ins = galaxy_instance(insfile, name = insname)
        gi = login(ins['key'],url=ins['url'])
        input_file_path = 'data_test/wf_inputs.yaml'
        params_file_path = 'data_test/wf_parameters.yaml'
        workflow = "MSA ClustalW"
        history_name = "history_test_input_files"
        invocation = run_workflow_files(gi,workflow,input_file_path,params_file_path,history_name)
        state = invocation_errors(gi,invocation)
        self.assertEqual(state,'ok')

    def test_galaxy_integration(self):
        insfile = 'data_test/parsec_creds.yaml'
        insname = 'local'
        fn = 'data_test/ls_orchid.fasta'
        resource_params = {
            'file' : insfile,
            'name' : insname
        }
        job = JobExecutorAtGalaxy()
        job.set_resource(resource_params)
        job.connect()
        # history_id = job.create_job_workspace(name = '22-10-2020')
        input_file_path = 'data_test/wf_inputs.yaml'
        params_file_path = 'data_test/wf_parameters.yaml'
        param_data = read_yaml_file(params_file_path)
        inputs_data = read_yaml_file(input_file_path)
        params = {
            'workflow' : "MSA ClustalW",
            'inputs': inputs_data,
            'parameters': param_data
        }
        history_id,invocation_id = job.submit('22-10-2020',params)
        self.assertIsNotNone(invocation_id,'error_at_invoke')
        status = job.job_status(invocation_id)
        print(status)
        job.remove_job_workspace(history_id)
        job.get_resuts(invocation_id)

    def test_wf_inst_2_inst(self):
        insfile = 'data_test/parsec_creds.yaml'
        instanceIn = instance(insfile,'beauvoir')
        instanceOut = instance(insfile,'balder')
        wf_in = get_workflow_from_name('MSA Clustaw')
        wf_dic_in = export_workflow(instanceIn, wf_in['id'])
        wf_out = import_workflow(instanceOut, wf_dic_in)
        wf_dic_out = export_workflow(wf_out)
        list_of_tools = check_tools(wf_dic_in,wf_dic_out)
        if isinstance(list_of_tools,list):
            tools_message = install_tools(instanceOut,list_of_tools)
            print(tools_message)

    def test_import_install(self):
        insfile = 'data_test/parsec_creds.yaml'
        insname = 'local'
        ins = galaxy_instance(insfile, name=insname)
        gi = login(ins['key'], url=ins['url'])
        workflow_path = '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/biobarcoding/workflows/Galaxy-Workflow-MSA_ClustalW.ga'
        # project_path = 'home/paula/Documentos/NEXTGENDEM/'
        wf = gi.workflows.import_workflow_from_local_path(workflow_path)
        import json
        with open(workflow_path, 'r') as f:
            wf_dict_in = json.load(f)
        wf_dict_out = gi.workflows.export_workflow_dict(wf['id'])
        list_of_tools = check_tools(wf_dict_out, wf_dict_in)
        install_tools(gi, list_of_tools)
        tools_message = install_tools(gi, list_of_tools)
        print(tools_message)

    def test_convert(self):
        import json
        input_path = '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/mrbayes_galaxy.json'
        with open(input_path, 'r') as f:
            galaxy_dict_in = json.load(f)
        inputs = galaxy_dict_in['inputs']
        convert = ToFormlyConverter()
        formly_json = convert.get_formly_json(inputs)
        json = '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/mrbayes_formly.json'
        with open(json, 'w') as file:
            file.write(formly_json)

    def test_convert_workflows(self):
        wfdict = {'clustalw': '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/clustalw_galaxy.json'}
        new_form_path = '/biobarcoding/inputs_schema/clustalw_wf_formly.json'
        convertToFormly(wfdict,new_form_path)
        path = '/biobarcoding/inputs_schema/clustalw_wf_formly.json'
        with open(path, 'w') as file:
            file.write(json)





if __name__ == '__main__':
    # MyTestCase.test_upload_file()
    # MyTestCase.test_run_workflow_from_file()
    # MyTestCase.test_change_parameter()
    # MyTestCase.test_inputs_files()
    # MyTestCase.test_import_install()
    MyTestCase.test_convert()
