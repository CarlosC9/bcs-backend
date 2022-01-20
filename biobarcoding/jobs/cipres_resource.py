import os
from shutil import rmtree
import requests
from python_cipres import client as cipres

from .ssh_process_adaptors import SSHProcessAdaptor
from .. import get_global_configuration_variable
from ..jobs import JobExecutorAtResource


class JobExecutorWithCipres(JobExecutorAtResource):

    def __init__(self, identity_job_id, create_local_workspace=True):
        super().__init__(identity_job_id, create_local_workspace)
        self.base_url = None
        self.username = None
        self.password = None
        self.app_name = None
        self.appID = None
        self.cipres_client = None
        self.last_job_remotely = True

        # RESOURCE

    def set_resource(self, resource_params):
        self.base_url = resource_params["jm_location"]["baseUrl"]
        self.username = resource_params["jm_credentials"]['username']
        self.password = resource_params["jm_credentials"]['password']
        self.app_name = resource_params["jm_credentials"]['appName']
        self.appID = resource_params["jm_credentials"]['appID']

    def check(self):
        return requests.head(self.base_url) == 200

    def connect(self):
        os.environ["URL"] = self.base_url
        os.environ["PASSWORD"] = self.password
        os.environ["KEY"] = self.appID
        self.cipres_client = cipres.Client(self.app_name, self.appID, self.username, self.password, self.base_url)

    def disconnect(self):
        pass

    # JOB EXECUTION
    def create_job_workspace(self):
        pass

    def job_workspace_exists(self):
        """In cipres there is no remote workspace per se. This function should return True, if the job
        hasn't been submitted or if the job has terminated and it has been deleted from CIPRES."""
        if os.path.exists(os.path.join(self.local_workspace, "job_handle.txt")):
            with open(os.path.join(self.local_workspace, "job_handle.txt", "r")) as f:
                job_handle = f.read()
                for job_status in self.cipres_client.listJobs(job_handle):
                    if job_status.jobHandle == job_handle:
                        return False

        return True

    def remove_job_workspace(self):
        """In CIPRES there is no remote workspace per se. This function delete the job records from CIPRES"""
        if os.path.exists(os.path.join(self.local_workspace, "job_handle.txt")):
            with open(os.path.join(self.local_workspace, "job_handle.txt", "r")) as f:
                job_handle = f.read()
            job_status = self.cipres_client.getjobStatus(job_handle)
            job_status.delete()

    def exists(self, job_context):
        return True

    def upload_file(self, job_context):
        pass

    def download_file(self, job_context):
        i = job_context["state_dict"]["idx"]
        if os.path.exists(os.path.join(self.local_workspace, "job_handle.txt")):
            with open(os.path.join(self.local_workspace, "job_handle.txt", "r")) as f:
                job_handle = f.read()
        else:
            return ""
        job_status = self.cipres_client.getjobStatus(job_handle)
        filename = self.get_download_files_list(job_context)[i]["file"]
        file_url = job_status.listResults()[filename].getUrl()
        cmd = (f"(nohup curl -u {self.username}:{self.password} -H cipres-appkey:{self.appID} " +
               f"-O -J --output-dir {self.local_workspace} {file_url}"
               f">>{self.log_filenames_dict['download_stdout']} </dev/null 2>>{self.log_filenames_dict['download_stderr']} " +
               f"& echo $!; wait $!; echo $? >> {self.local_workspace}/$!.exit_status)")

        print(cmd)
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().rstrip()
        print(f"PID: {pid}")
        return pid

    def submit(self, process):
        params = process["inputs"]["vParams"]
        input_files = process["inputs"]["inputParams"]
        metadata = {}
        try:
            job_status = self.cipres_client.submitJob(params, input_files, metadata, validateOnly=True)  # TODO
            with open(os.path.join(self.local_workspace, "job_handle.txt"), "w") as f:
                f.write(job_status.jobHandle)
            print(f"JobHandle {job_status.jobHandle}")
            return job_status.jobHandle
        except Exception as e:
            print(e)
            return ""  # error

    def step_status(self, job_context):
        state = job_context["state_dict"]["state"]
        if state in ["prepare", "transfer_data"]:
            status = "ok"
        else:
            pid = job_context.get("pid")  # pid is jobHandle in case previous_state = submit
            if pid is None:
                status = "none"
            elif pid == "":
                status = ""  # error
            elif state == "submit":
                job_status = self.cipres_client.getjobStatus(pid)
                if job_status.isError:
                    status = ""
                elif job_status.isDone:
                    status = "ok"
                else:
                    status = "running"
                if status == "ok" or status == "":
                    self.write_remote_logs(job_context["state_dict"], job_status)
            else:
                status = self.local_job_status(pid)

        return status

    def get_cipres_logs(self, job_status):
        if not job_status.jobHandle and job_status.commandline:
            return "Submission validated.  Commandline is: '%s'" % job_status.commandline

        s = "Job=%s" % job_status.jobHandle
        if job_status.terminalStage:
            if job_status.failed:
                s += ", failed at stage %s" % job_status.jobStage
            else:
                s += ", finished, results are at %s" % job_status.resultsUrl
        else:
            s += ", not finished, stage=%s" % job_status.jobStage
        s += "\n MESSAGES \n"
        for m in job_status.messages:
            s += "\t%s\n" % m
        if job_status.metadata:
            s += "METADATA\n"
            for key in job_status.metadata:
                s += "\t%s=%s\n" % (key, job_status.metadata[key])
        else:
            s += "There is no metadata\n"

        return s

    def write_remote_logs(self, state_dict, job_status):
        if state_dict["state"] == "submit":
            logs = self.get_cipres_logs(job_status)
            if job_status.isDone:
                with open(self.log_filenames_dict['submit_stdout'], "a") as file:
                    file.write(logs)
            else:
                with open(self.log_filenames_dict['submit_stderr'], "a") as file:
                    file.write(logs)

    def cancel_job(self, native_id):
        #native_id is the job handle in CIPRES
        job_status = self.cipres_client.getjobStatus(native_id)
        job_status.delete()

    def get_upload_files_list(self, job_context):
        return []

    def get_download_files_list(self, job_context):
        return job_context["results"]

    def check_resource(self):
        """
        Can connect
        CPU and GPU use
        Available storage space

        :return:
        """
        # TODO:
        pass
