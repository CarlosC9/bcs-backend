import os
import abc

from biobarcoding.jobs.process_adaptor import ProcessAdaptor

SCRIPT_KEY = "script"
SCRIPT_FILES_KEY = "script_files"
SCRIPT_PARAMS_KEY = "script_params"
RESULTS_FILES_KEY = "result_files"


class SSHProcessAdaptor(ProcessAdaptor, abc.ABC):
    '''The methods of the subinterfaces are all private'''

    ASSETS_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "ssh_process_assets")

    def adapt_job_context(self, job_context):
        process_parameters = job_context["process"]["inputs"]["parameters"]
        print(process_parameters)
        new_process_parameters = {
            SCRIPT_KEY: self._get_script_filename(),
            SCRIPT_FILES_KEY: self._get_script_files_list() +
                              [{
                                  "remote_name": self._get_script_filename(),
                                  "file": os.path.join(self.ASSETS_FOLDER, self._get_script_filename()),
                                  "type": "sh"
                              }],
            SCRIPT_PARAMS_KEY: self._get_script_params_string(process_parameters),
            RESULTS_FILES_KEY: self._get_results_files_list(process_parameters),
        }
        job_context["process"]["inputs"]["parameters"] = new_process_parameters
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
    def _get_results_files_list(self, process_parameters):
        raise NotImplementedError
