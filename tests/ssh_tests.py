import asyncio
import json
import os
import unittest

from biobarcoding.jobs.ssh_resource import JobExecutorWithSSH


class SSHTestCase(unittest.TestCase):

    def upload_files(self, remote_client):
        remote_client.upload_file("", "/home/daniel/Documentos/GIT/bcs-backend/tests/ssh_data_test/add_line.py", "add_line.py")
        remote_client.upload_file("", "/home/daniel/Documentos/GIT/bcs-backend/tests/ssh_data_test/count_lines.py", "count_lines.py")
        remote_client.upload_file("", "/home/daniel/Documentos/GIT/bcs-backend/tests/ssh_data_test/myfile.txt", "myfile.txt")
        remote_client.upload_file("", "/home/daniel/Documentos/GIT/bcs-backend/tests/ssh_data_test/test.sh", "test.sh")

    def main_test(self):
        ssh_data_dirname = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ssh_data_test")
        with open(os.path.join(ssh_data_dirname, 'test_params.json')) as json_file:
            params = json.load(json_file)

        remote_client = JobExecutorWithSSH()
        remote_client.set_resource(params)
        remote_client.connect()
        remote_client.create_job_workspace(params['remote_workspace'])
        self.upload_files(remote_client)
        pid = remote_client.submit(params['remote_workspace'], params)
        print(remote_client.job_status(pid))
        remote_client.create_job_workspace("results")
        remote_client.move_file("myfile.txt", "results/myfile.txt")
        remote_client.retrieve_file("results/myfile.txt", "/home/daniel/Documentos")
        remote_client.retrieve_directory("results/", "/home/daniel/Documentos")
        remote_client.remove_job_workspace(params['remote_workspace'])
        remote_client.disconnect()

    def killed_test(self):
        remote_client = JobExecutorWithSSH()
        ssh_data_dirname = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ssh_data_test")
        with open(os.path.join(ssh_data_dirname, 'test_killed_params.json')) as json_file:
            params = json.load(json_file)
            remote_client.set_resource(params)

        if remote_client.check():
            remote_client.connect()
            remote_client.make_directory(params["remote_workspace"])
            remote_client.upload_file("", os.path.join(ssh_data_dirname, "test_killed.sh"), "test_killed.sh")
            if remote_client.check_files(["test_killed.sh"]):
                asyncio.get_event_loop().run_until_complete(remote_client.submit(params['remote_workspace'], params))
                status = remote_client.job_status("")
                print(f"Job status running: {status}")
                remote_client.cancel_job("")
                print(f"Job Killed")
                remote_client.job_status("")
                remote_client.remove_file("test_killed.sh")
                file_exists = remote_client.check_files(["test_killed"])
                print(f"File doesn't exist: {not file_exists}")
                remote_client.remove_job_workspace(params["remote_workspace"])


if __name__ == '__main__':
    SSHTestCase().main_test()
