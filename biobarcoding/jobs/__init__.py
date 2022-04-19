import abc
import os
import pathlib
from abc import ABC
from typing import Dict

from ..tasks import celery_app
from .. import get_global_configuration_variable


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
    LOCAL_WORKSPACE = os.sep + "tmp"  # TODO read from config file
    prefix = "jobs"
    LOG_FILENAMES_DICT = {
        "prepare_stdout": f"{prefix}.prepare.stdout.log",
        "prepare_stderr": f"{prefix}.prepare.stderr.log",
        "export_stdout": f"{prefix}.export.stdout.log",
        "export_stderr": f"{prefix}.export.stderr.log",
        "upload_stdout": f"{prefix}.upload.stdout.log",
        "upload_stderr": f"{prefix}.upload.stderr.log",
        "submit_stdout": f"{prefix}.submit.stdout.log",
        "submit_stderr": f"{prefix}.submit.stderr.log",
        "download_stdout": f"{prefix}.download.stdout.log",
        "download_stderr": f"{prefix}.download.stderr.log",
        "store_stdout": f"{prefix}.store.stdout.log",
        "store_stderr": f"{prefix}.store.stderr.log",
        "cleanup_stdout": f"{prefix}.cleanup.stdout.log",
        "cleanup_stderr": f"{prefix}.cleanup.stderr.log",
        "universal_log": f"{prefix}.stdout.log",
    }

    def __init__(self, identity_job_id, create_local_workspace=True):
        self.local_workspace = os.path.join(
            get_global_configuration_variable("JOBS_LOCAL_WORKSPACE", self.LOCAL_WORKSPACE),
            identity_job_id)
        if create_local_workspace:
            if not os.path.exists(self.local_workspace):
                pathlib.Path(self.local_workspace).mkdir(parents=True, exist_ok=True)

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

    def get(self, job_context, create_local_workspace=True):
        job_executor_name = job_context["resource"].get("jm_type")
        resource_param = job_context["resource"]
        k = (job_executor_name, resource_param["name"])
        job_id = job_context.get("job_id")
        if k not in self.execs:
            self.execs[k] = JobExecutorAtResourceFactory._create(job_executor_name,
                                                                 resource_param,
                                                                 create_local_workspace=create_local_workspace,
                                                                 job_id=job_id,
                                                                 identity_id=job_context.get("identity_id"))
        return self.execs[k]

    @staticmethod
    def _create(job_executor_name: str, resource_param: Dict, create_local_workspace, **kwargs):
        if job_executor_name.lower() == "galaxy":
            from biobarcoding.jobs.galaxy_resource import JobExecutorAtGalaxy
            tmp = JobExecutorAtGalaxy(str(kwargs["job_id"]) + "_" + str(kwargs["identity_id"]))
            tmp.set_resource(resource_param)
            return tmp
        elif job_executor_name.lower() == "ssh":
            from ..jobs.ssh_resource import JobExecutorWithSSH
            tmp = JobExecutorWithSSH(str(kwargs["job_id"]) + "_" + str(kwargs["identity_id"]), create_local_workspace)
            tmp.set_resource(resource_param)
            tmp.connect()
            return tmp
        elif job_executor_name.lower() == "slurm":
            from .slurm_ssh_resource import JobExecutorWithSlurm
            remote_workspace = None
            if resource_param["name"].lower().startswith("teide"):
                remote_workspace = get_global_configuration_variable("TEIDE_REMOTE_WORKSPACE")
            tmp = JobExecutorWithSlurm(str(kwargs["job_id"]) + "_" + str(kwargs["identity_id"]),
                                       create_local_workspace,
                                       remote_workspace)
            tmp.set_resource(resource_param)
            tmp.connect()
            return tmp
        elif job_executor_name.lower() == "cipres":
            from .cipres_resource import JobExecutorWithCipres
            tmp = JobExecutorWithCipres(str(kwargs["job_id"]) + "_" + str(kwargs["identity_id"]), create_local_workspace)
            tmp.set_resource(resource_param)
            tmp.connect()
            return tmp

