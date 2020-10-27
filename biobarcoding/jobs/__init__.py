from abc import ABC
from celery import chain
from typing import Dict

from biobarcoding.tasks.definitions import *


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
    def set_resource(self, params):
        pass

    def check(self):
        pass

    def connect(self):
        pass

    def disconnect(self):
        pass

    # JOB EXECUTION
    def set_credentials(self, credentials):
        """ Different from connecting to the resource, job submission may require identifying user, to check
            priority and the like
            """
        pass

    def get_quotas_for_current_credentials(self):
        pass

    def create_job_workspace(self, name):
        pass

    def remove_job_workspace(self, name):  # After Job is completed (or if Job was not started)
        pass

    def upload_file(self, workspace, local_filename, remote_location):
        pass

    def move_file(self, remote_source, remote_destination):
        pass

    def remove_file(self, remote_filename):
        pass

    def submit(self, workspace, params):
        pass

    def job_status(self, native_id):
        pass

    def cancel_job(self, native_id):
        pass


class JobExecutorAtResourceFactory:
    def __init__(self):
        self.execs = dict()

    def get(self, job_executor_name, resource_param: Dict):
        k = (job_executor_name, resource_param["name"])
        if k not in self.execs:
            self.execs[k] = JobExecutorAtResourceFactory._create(job_executor_name, resource_param)

        return self.execs[k]

    @staticmethod
    def _create(job_executor_name: str, resource_param: Dict):
        if job_executor_name.lower() == "galaxy":
            from biobarcoding.jobs.galaxy_resource import JobExecutorAtGalaxy
            tmp = JobExecutorAtGalaxy()
            tmp.set_resource(resource_param["jm_location"].update(resource_param["jm_credentials"]))
            return tmp