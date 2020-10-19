from abc import ABC
from celery import chain
from biobarcoding.tasks.definitions import *


class JobManagementAPI:
    @staticmethod
    def submit(params):
        """ In the params there must be a reference to the resource.
            Also the desired job executor (a resource can have multiple -although it shouldn't-)
            Check if there is job executor for the resource.
            If there exists, use it. If not, create one

         """
        # TODO Which process?
        #
        params2 = params
        chain(wf1_prepare_workspace.s(params2),
              wf1_export_to_supported_file_formats.s(),
              wf1_transfer_data_from_resource.s(),
              wf1_submit.s(),
              wf1_wait_until_execution_starts.s(),
              wf1_wait_for_execution_finish.s(),
              wf1_transfer_data_from_resource.s(),
              wf1_import_into_database.s(),
              wf1_cleanup_workspace.s(),
              wf1_complete_succesfully.s()).delay()

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

    def job_status(self, id):
        pass

    def cancel_job(self, id):
        pass
