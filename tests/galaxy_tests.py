import unittest
from pathlib import Path
from bioblend import galaxy
from biobarcoding.jobs.galaxy_resource import *

class MyTestCase(unittest.TestCase):
    def test_upload_file(self):
        user_key = 'af107bf81f146b6746944b9488986822'
        url = 'http://127.0.0.1:8080'
        gi = galaxy.GalaxyInstance(url = url, key = user_key)
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
        user_key = 'af107bf81f146b6746944b9488986822'
        url = 'http://127.0.0.1:8080'
        gi = login(user_key,url=url)
        fn = 'data_test/ls_orchid.fasta'
        invocation = run_workflow(gi,"Workflow_Input",fn,'marK1')
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
        user_key = 'af107bf81f146b6746944b9488986822'
        url = 'http://127.0.0.1:8080'
        gi = login(user_key, url=url)
        fn = 'data_test/ls_orchid.fasta'
        workflow = "Workflow_Input"
        step = '2'
        parameter_name = 'equi_freq'
        value = '"e"'
        w_id = workflow_id(gi,workflow)
        new_parameters = set_parameters(gi,w_id,step,parameter_name,value) # esta función me la podría cargar
        invocation = run_workflow(gi, workflow, fn, 'marK1', params=new_parameters)
        state = invocation['state']

        self.assertEqual(state, 'new', 'invocation failed')
        while invocation_errors(gi, invocation) == 0 and invocation_percent_complete(gi, invocation) < 100:
            pass
        errors = invocation_errors(gi, invocation)
        step_job = get_job(gi,invocation,step)
        self.assertEqual(errors, 0, 'There is an error in the invocation')
        self.assertEqual(step_job['params'][parameter_name].strip('"'),new_parameters[step][parameter_name])
        # TODO it seems that it worked by looking in galaxy GUI parameter is still ML


if __name__ == '__main__':
    # MyTestCase.test_upload_file()
    # MyTestCase.test_run_workflow_from_file()
    MyTestCase.test_change_parameter()
