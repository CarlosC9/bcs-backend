import os
import requests
import subprocess
from python_cipres import client as cipres
from datetime import datetime

from .cipres_process_adaptors import CipresProcessAdaptor
from ..jobs import JobExecutorAtResource

######## OVERRIDE CIPRES METHOD FOR GETTING THE STAGE OF THE MESSAGES #######

class MyJobStatus(cipres.JobStatus):

    def __init__(self, client, jobUrl=None, xml=None):
        super(MyJobStatus, self).__init__(client, jobUrl, xml)
        # Override note: self.myJobStage better than self.jobStage
        self.myJobStage = None

    def __parseJobStatus__(self, xml):
        if xml.find("commandline") is not None:
            self.commandline = xml.find("commandline").text
        if xml.find("selfUri") is not None:
            self.jobUrl = xml.find("selfUri").find("url").text
        if xml.find("jobHandle") is not None:
            self.jobHandle = xml.find("jobHandle").text
        if xml.find("jobStage") is not None:
            self.jobStage = xml.find("jobStage").text
        if xml.find("terminalStage") is not None:
            self.terminalStage = (xml.find("terminalStage").text == "true")
        if xml.find("failed") is not None:
            self.failed = (xml.find("failed").text == "true")
        if xml.find("resultsUri") is not None:
            self.resultsUrl = xml.find("resultsUri").find("url").text
        if xml.find("workingDirUri") is not None:
            self.workingDirUrl = xml.find("workingDirUri").find("url").text
        if xml.find("dateSubmitted") is not None:
            self.dateSubmitted = xml.find("dateSubmitted").text
        if xml.find("messages") is not None:
            # Override note: adding final stage in message to self.myJobStage
            last_date = None
            for m in xml.find("messages"):
                date = datetime.strptime(m.find("timestamp").text.rsplit("-", 1)[0], "%Y-%m-%dT%H:%M:%S")
                if last_date is None or date > last_date:
                    last_date = date
                    self.myJobStage = m.find("stage").text
                    print("-----------------TESTING JOB STAGE-----------------")
                    print(m)
                    print(self.myJobStage)
                self.messages.append("%s: %s" % (m.find("timestamp").text, m.find("text").text))
        if xml.find("metadata") is not None:
            for e in xml.find("metadata").findall("entry"):
                self.metadata[e.find("key").text] = e.find("value").text


# Override cipres.JobStatus to have self.myJobStatus
cipres.JobStatus = MyJobStatus

####### END OF OVERRIDE #######


class JobExecutorWithCipres(JobExecutorAtResource):

    JOB_STATES_DICT = {
        'QUEUE': "running",
        'COMMANDRENDERING': "running",
        'INPUTSTAGING': "running",
        'SUBMITTED': "running",
        'LOAD_RESULTS': "running",
        'COMPLETED': "ok"
    }

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
        return requests.head(self.base_url).status_code == 200

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
            with open(os.path.join(self.local_workspace, "job_handle.txt"), "r") as f:
                job_handle = f.read()
            cmd = f"curl -s -i -u {self.username}:{self.password} -H cipres-appkey:{self.appID} " +\
                  f"{self.base_url}/job/{self.username}/{job_handle} | grep 'Job not found'"
            result = subprocess.run(cmd, shell=True)
            print(cmd)
            print(result.returncode)
            return result.returncode == 1

        return True

    def remove_job_workspace(self):
        """In CIPRES there is no remote workspace per se. This function delete the job records from CIPRES"""
        if os.path.exists(os.path.join(self.local_workspace, "job_handle.txt")):
            with open(os.path.join(self.local_workspace, "job_handle.txt"), "r") as f:
                job_handle = f.read()
            cmd = f"curl -u {self.username}:{self.password} -H cipres-appkey:{self.appID} " + \
                  f"-X DELETE {self.base_url}/job/{self.username}/{job_handle}"
            r = subprocess.run(cmd, shell=True)
            print(cmd)
            print(r)

    def exists(self, job_context):
        i = job_context["state_dict"]["idx"]
        if job_context["state_dict"]["state"] == "upload":
            filename = self.get_upload_files_list(job_context)[i].get("remote_name")
        else:
            filename = self.get_download_files_list(job_context)[i].get("file")
        path = os.path.join(self.local_workspace, filename)
        return os.path.exists(path)

    def upload_file(self, job_context):
        '''The upload actually is a copy of the script files to the local workspace'''
        i = job_context["state_dict"]["idx"]
        path = self.get_upload_files_list(job_context)[i].get("file")
        target_path = os.path.join(self.local_workspace, self.get_upload_files_list(job_context)[i].get("remote_name"))
        cmd = (f"(nohup bash -c \"cp {path} {target_path}\" >>{self.log_filenames_dict['upload_stdout']} " +
               f"</dev/null 2>>{self.log_filenames_dict['upload_stderr']} & echo $!; wait $!; echo $? >> " +
               f"{self.local_workspace}/$!.exit_status)")
        print(cmd)
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().rstrip()
        print(f"PID: {pid}")
        return pid

    def download_file(self, job_context):
        i = job_context["state_dict"]["idx"]
        if os.path.exists(os.path.join(self.local_workspace, "job_handle.txt")):
            with open(os.path.join(self.local_workspace, "job_handle.txt"), "r") as f:
                job_handle = f.read()
        else:
            return ""
        job_status = self.cipres_client.getJobStatus(job_handle)
        print("------------------------TESTING DOWNLOAD----------------------------")
        print(job_status.listResults())
        filename = self.get_download_files_list(job_context)[i]["file"]
        remote_filename = self.get_download_files_list(job_context)[i]["remote_name"]
        file_url = job_status.listResults()[remote_filename].getUrl()
        cmd = (f"(nohup curl -o {self.local_workspace}/{filename} -u {self.username}:{self.password} " +
               f"-H cipres-appkey:{self.appID} -O -J {file_url}"
               f">>{self.log_filenames_dict['download_stdout']} </dev/null 2>>{self.log_filenames_dict['download_stderr']} " +
               f"& echo $!; wait $!; echo $? >> {self.local_workspace}/$!.exit_status)")

        print(cmd)
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().rstrip()
        print(f"PID: {pid}")
        return pid

    def submit(self, process):
        params = process["inputs"]["adapted_parameters"]
        script_file = params[CipresProcessAdaptor.SCRIPT_KEY][0]["remote_name"]
        script_params = params[CipresProcessAdaptor.SCRIPT_PARAMS_KEY]
        credentials = f"username={self.username} base_url={self.base_url} password={self.password} appID={self.appID} app_name={self.app_name} "
        print("--------------------------------- SUBMIT ------------------------------------")
        print(script_params)
        env_variables = credentials + script_params
        print("------------------------------------------------------")
        print(env_variables)
        cmd = (
                f"cd {self.local_workspace} && chmod +x {script_file} &> /dev/null " +
                f"&& export {env_variables} && (cd {self.local_workspace} && nohup ./{script_file} " +
                f">/{self.local_workspace}/{os.path.basename(self.LOG_FILENAMES_DICT['submit_stdout'])} " +
                f"</dev/null 2>/{self.local_workspace}/{os.path.basename(self.LOG_FILENAMES_DICT['submit_stderr'])}" +
                f"& echo $!; wait $!; echo $? >> {self.local_workspace}/$!.exit_status)")

        print(repr(cmd))
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().rstrip()
        self.last_job_remotely = True
        print(f"PID: {pid}")
        return pid

    def step_status(self, job_context):
        state = job_context["state_dict"]["state"]
        if state in ["prepare"]:
            status = "ok"
        else:
            pid = job_context.get("pid")
            if pid is None:
                status = "none"
            elif pid == "":
                status = ""  # error
            elif state == "submit":
                status = self.local_job_status(pid)
                if status == "ok":
                    if os.path.exists(os.path.join(self.local_workspace, "job_handle.txt")):
                        with open(os.path.join(self.local_workspace, "job_handle.txt"), "r") as f:
                            job_handle = f.read()
                        job_status = self.cipres_client.getJobStatus(job_handle)
                        '''if job_status.isError:
                            print("aqui")
                            print(job_status.myJobStage)
                            print(job_status.jobStage)
                            status = ""
                        elif job_status.isDone:
                            print("aqui 2")
                            status = "ok"
                        else:'''
                        print("aqui 3")
                        print(job_status.myJobStage)
                        print(job_status.jobStage)
                        status = self.JOB_STATES_DICT[job_status.jobStage]
                        if status == "ok" or status == "":
                            self.write_remote_logs(job_context["state_dict"], job_status)
                    else:
                        status = ""
            else:
                status = self.local_job_status(pid)

        return status

    def get_cipres_logs(self, job_status):
        if not job_status.jobHandle and job_status.commandline:
            return "Submission validated.  Commandline is: '%s'" % job_status.commandline

        s = "Job=%s" % job_status.jobHandle
        if job_status.terminalStage:
            if job_status.failed:
                s += ", failed at stage %s" % job_status.myJobStage
            else:
                s += ", finished, results are at %s" % job_status.resultsUrl
        else:
            s += ", not finished, stage=%s" % job_status.myJobStage
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
            logs = "\n" + self.get_cipres_logs(job_status)
            if job_status.isDone:
                with open(self.log_filenames_dict['submit_stdout'], "a") as file:
                    file.write(logs)
            else:
                with open(self.log_filenames_dict['submit_stderr'], "a") as file:
                    file.write(logs)

    def cancel_job(self, native_id):
        os.kill(native_id)
        if os.path.exists(os.path.join(self.local_workspace, "job_handle.txt")):
            with open(os.path.join(self.local_workspace, "job_handle.txt"), "r") as f:
                job_handle = f.read()
            job_status = self.cipres_client.getJobStatus(job_handle)
            job_status.delete()

    def get_upload_files_list(self, job_context):
        return list(job_context["process"]["inputs"]["adapted_parameters"][CipresProcessAdaptor.SCRIPT_FILES_KEY]) + \
               list(job_context["process"]["inputs"]["adapted_parameters"].get("scripts", []))

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
