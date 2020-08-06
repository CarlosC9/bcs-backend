import unittest
from pathlib import Path
from bioblend import galaxy
from biobarcoding.jobs.galaxy_resource import login,run_workflow,download_result,list_invocation_results

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
        # TODO control job execution
        #results = list_invocation_results(gi,invocation_id=invocation['id'])
        #download_result(gi,results,'data_test/')
        #self.assertIsInstance(results,list,'There are no results')



if __name__ == '__main__':
    # MyTestCase.test_upload_file()
    MyTestCase.test_run_workflow_from_file()
