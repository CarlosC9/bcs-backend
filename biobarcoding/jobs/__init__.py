from abc import ABC


class JobManagementAPI:
    def submit(self, params):
        """ In the params there must be a reference to the resource.
            Also the desired job executor (a resource can have multiple -although it shouldn't-)
            Check if there is job executor for the resource.
            If there exists, use it. If not, create one

         """

    def check(self, job_id):
        pass

    def cancel(self, job_id):
        pass

    def list(self, filter):
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
