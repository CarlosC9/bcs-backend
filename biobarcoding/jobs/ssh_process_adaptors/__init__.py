import os
from abc import ABC

from biobarcoding.jobs.process_adaptor import ProcessAdaptor


class SSHProcessAdaptor(ProcessAdaptor, ABC):
    '''The methods of the subinterfaces are all private'''

    ASSETS_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "ssh_process_adaptors",
                                 "ssh_process_assets")
    SCRIPT_KEY = "script"
    SCRIPT_FILES_KEY = "script_files"
    SCRIPT_PARAMS_KEY = "script_params"
    RESULTS_FILES_KEY = "result_files"

    def adapt_job_context(self, job_context):
        process_parameters = job_context["process_params"]["parameters"]
        new_process_parameters = {
            self.SCRIPT_KEY: self._get_script_filename(),
            self.SCRIPT_FILES_KEY: self._get_script_files_list(),
            self.SCRIPT_PARAMS_KEY: self._get_script_params_string(process_parameters),
            self.RESULTS_FILES_KEY: self._get_results_files_list(),
        }
        job_context["process_params"]["parameters"] = new_process_parameters
        return job_context

    @abc.abstractmethod
    def _get_script_filename(self):
        raise NotImplementedError

    @abc.abstractmethod
    def _get_script_files_list(self):
        raise NotImplementedError

    @abc.abstractmethod
    def _get_script_params_string(self, process_parameters):
        raise NotImplementedError

    @abc.abstractmethod
    def _get_results_files_list(self):
        raise NotImplementedError