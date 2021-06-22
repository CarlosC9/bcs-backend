import os
from abc import ABC
import abc
from typing import Dict

from biobarcoding.rest.file_manager import FilesAPI
from biobarcoding.tasks import celery_app


class JobManagementAPI:
    @staticmethod
    def submit(job_context):
        """ In the params there must be a reference to the resource.
            Also the desired job executor (a resource can have multiple -although it shouldn't-)
            Check if there is job executor for the resource.
            If there exists, use it. If not, create one

         """
        celery_app.signature("prepare").delay(job_context)

    @staticmethod
    def check(job_id):
        """
        Check the status of a job

        :param job_id:
        :return:
        """
        pass

    @staticmethod
    def cancel(job_id):
        pass

    @staticmethod
    def list(filter):
        pass

    @staticmethod
    def get(job_id):
        pass


class JobExecutorAtResource(ABC):
    # TODO acordar el path con Rafa
    LOCAL_WORKSPACE = os.sep + "tmp"# TODO leer desde config file
    LOG_FILENAMES_DICT = {
        "prepare_stdout": "bcs.prepare.stdout.log",
        "prepare_stderr": "bcs.prepare.stderr.log",
        "export_stdout": "bcs.export.stdout.log",
        "export_stderr": "bcs.export.stderr.log",
        "upload_stdout": "bcs.upload.stdout.log",
        "upload_stderr": "bcs.upload.stderr.log",
        "submit_stdout": "bcs.submit.stdout.log",
        "submit_stderr": "bcs.submit.stderr.log",
        "download_stdout": "bcs.download.stdout.log",
        "download_stderr": "bcs.download.stderr.log",
        "store_stdout": "bcs.store.stdout.log",
        "store_stderr": "bcs.store.stderr.log",
        "cleanup_stdout": "bcs.cleanup.stdout.log",
        "cleanup_stderr": "bcs.cleanup.stderr.log",
        "universal_log": "bcs.stdout.log",
    }

    def __init__(self, identity_job_id):
        self.local_workspace = os.path.join(self.LOCAL_WORKSPACE, identity_job_id)
        if not os.path.exists(self.local_workspace):
            os.mkdir(self.local_workspace)
        self.log_filenames_dict = {}
        for k, v in self.LOG_FILENAMES_DICT.items():
            self.log_filenames_dict[k] = os.path.join(self.local_workspace, v)

    # RESOURCE
    @abc.abstractmethod
    def set_resource(self, params):
        raise NotImplementedError

    @abc.abstractmethod
    def check(self):
        raise NotImplementedError

    @abc.abstractmethod
    def connect(self):
        raise NotImplementedError

    @abc.abstractmethod
    def disconnect(self):
        raise NotImplementedError

    # JOB EXECUTION
    @abc.abstractmethod
    def create_job_workspace(self):
        raise NotImplementedError

    @abc.abstractmethod
    def job_workspace_exists(self):
        raise NotImplementedError

    @abc.abstractmethod
    def remove_job_workspace(self):  # After Job is completed (or if Job was not started)
        raise NotImplementedError

    @abc.abstractmethod
    def exists(self, job_context):
        raise NotImplementedError

    @abc.abstractmethod
    def upload_file(self, job_context):
        raise NotImplementedError

    @abc.abstractmethod
    def download_file(self, job_context):
        raise NotImplementedError

    @abc.abstractmethod
    def submit(self, job_context):
        raise NotImplementedError

    @abc.abstractmethod
    def step_status(self, job_context):
        raise NotImplementedError

    @abc.abstractmethod
    def cancel_job(self, native_id):
        raise NotImplementedError

    @abc.abstractmethod
    def get_upload_files_list(self, job_context):
        raise NotImplementedError

    @abc.abstractmethod
    def get_download_files_list(self, job_context):
        raise NotImplementedError

    def local_job_status(self, pid):
        exit_status = "none"
        if pid:
            if os.path.exists(f"{self.local_workspace}/{pid}.exit_status"):
                with open(f"{self.local_workspace}/{pid}.exit_status", "r") as f:
                    exit_status = f.readline().strip()
                    if exit_status.strip() == "0":
                        exit_status = "ok"
                    else:
                        print(f"Error executing get with pid: {pid}. Exit status = {exit_status}")
                        exit_status = ""  # This means error
            else:
                exit_status = "running"

        return exit_status


class JobExecutorAtResourceFactory:
    def __init__(self):
        self.execs = dict()

    def get(self, job_context):
        job_executor_name = job_context["resource"].get("jm_type")
        resource_param = job_context["resource"]
        k = (job_executor_name, resource_param["name"])
        if k not in self.execs:
            self.execs[k] = JobExecutorAtResourceFactory._create(job_executor_name,
                                                                 resource_param,
                                                                 job_id=job_context.get("job_id"),
                                                                 identity_id=job_context.get("identity_id"))
        return self.execs[k]

    @staticmethod
    def _create(job_executor_name: str, resource_param: Dict, **kwargs):
        if job_executor_name.lower() == "galaxy":
            from biobarcoding.jobs.galaxy_resource import JobExecutorAtGalaxy
            tmp = JobExecutorAtGalaxy(str(kwargs["job_id"]) + "_" + str(kwargs["identity_id"]))
            tmp.set_resource(resource_param)
            return tmp
        elif job_executor_name.lower() == "ssh":
            from biobarcoding.jobs.ssh_resource import JobExecutorWithSSH
            tmp = JobExecutorWithSSH(str(kwargs["job_id"]) + "_" + str(kwargs["identity_id"]))
            tmp.set_resource(resource_param)
            tmp.connect()
            return tmp
