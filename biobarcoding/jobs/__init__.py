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
    def create_job_workspace(self, name):
        raise NotImplementedError

    @abc.abstractmethod
    def remove_job_workspace(self, name):  # After Job is completed (or if Job was not started)
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
    def job_status(self, job_context):
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

    @abc.abstractmethod
    def get_store_path(self, job_context):
        raise NotImplementedError


class JobExecutorAtResourceFactory:
    def __init__(self):
        self.execs = dict()

    def get(self, job_context):
        job_executor_name = job_context["resource"].get("jm_type")
        resource_param = job_context["resource"]
        k = (job_executor_name, resource_param["name"])
        job_id = job_context.get("job_id")
        if k not in self.execs:
            self.execs[k] = JobExecutorAtResourceFactory._create(job_executor_name, resource_param, job_id=job_id)
        return self.execs[k]

    @staticmethod
    def _create(job_executor_name: str, resource_param: Dict, **kwargs):
        if job_executor_name.lower() == "galaxy":
            from biobarcoding.jobs.galaxy_resource import JobExecutorAtGalaxy
            tmp = JobExecutorAtGalaxy()
            tmp.set_resource(resource_param)
            return tmp
        elif job_executor_name.lower() == "ssh":
            from biobarcoding.jobs.ssh_resource import JobExecutorWithSSH
            tmp = JobExecutorWithSSH(str(kwargs["job_id"]))
            tmp.set_resource(resource_param)
            tmp.connect()
            return tmp
