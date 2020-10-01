import unittest
from pathlib import Path
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

    def test_run_workflow_from_file(self):
        insfile = 'data_test/parsec_creds.yaml'
        insname = 'beauvoir3'
        ins = galaxy_instance(insfile, name = insname)
        gi = login(ins['key'],url=ins['url'])
        fn = 'data_test/ls_orchid.fasta'
        workflow = 'MSA ClustalW'
        invocation = run_workflow(gi,workflow,fn,'marK1')
        state = invocation['state']
        self.assertEqual(state,'new','invocation failed')
        while invocation_errors(gi, invocation) == 0 and invocation_percent_complete(gi, invocation) < 100:
            pass
        errors = invocation_errors(gi,invocation)
        self.assertEqual(errors,0,'There is an error in the invocation')
        results = list_invocation_results(gi,invocation_id=invocation['id'])
        download_result(gi,results,'data_test/')
        self.assertIsInstance(results,list,'There are no results')

    def test_change_parameter(self):
        insfile = 'data_test/parsec_creds.yaml'
        insname = 'beauvoir3'
        ins = galaxy_instance(insfile, name = insname)
        gi = login(ins['key'],url=ins['url'])
        fn = 'data_test/ls_orchid.fasta'
        workflow = "PhyML_test_labels"
        step = '2'
        parameter_name = 'nbSubstCat'
        value = '3'
        w_id = workflow_id(gi,workflow)
        new_parameters = set_parameters(step,parameter_name,value)
        invocation = run_workflow(gi, workflow, fn, 'marK1', params=new_parameters)
        state = invocation['state']

        self.assertEqual(state, 'new', 'invocation failed')
        while invocation_errors(gi, invocation) == 0 and invocation_percent_complete(gi, invocation) < 100:
            pass
        errors = invocation_errors(gi, invocation)
        step_job = get_job(gi,invocation,step)
        self.assertEqual(errors, 0, 'There is an error in the invocation')
        self.assertEqual(step_job['params'][parameter_name].strip('"'),new_parameters[step][parameter_name])


    def test_inputs_files(self):
        insfile = 'data_test/parsec_creds.yaml'
        insname = 'beauvoir3'
        ins = galaxy_instance(insfile, name = insname)
        gi = login(ins['key'],url=ins['url'])
        input_file_path = 'data_test/wf_inputs.yaml'
        params_file_path = 'data_test/wf_parameters.yaml'
        workflow = "PhyML_test_labels"
        history_name = "history_test_input_files"
        invocation = run_workflow_files(gi,workflow,input_file_path,params_file_path,history_name)
        state = invocation['state']
        self.assertEqual(state, 'new', 'invocation failed')
        while invocation_errors(gi, invocation) == 0 and invocation_percent_complete(gi, invocation) < 100:
            pass
        errors = invocation_errors(gi, invocation)
        self.assertEqual(errors, 0, 'There is an error in the invocation')






if __name__ == '__main__':
    MyTestCase.test_upload_file()
    MyTestCase.test_run_workflow_from_file()
    MyTestCase.test_change_parameter()
    MyTestCase.test_inputs_files()
