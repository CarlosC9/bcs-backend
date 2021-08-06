import abc
import os

from ..process_adaptor import ProcessAdaptor

SCRIPT_KEY = "scripts"
SCRIPT_FILES_KEY = "scripts_files"
SCRIPT_PARAMS_KEY = "script_params"
RESULTS_KEY = "results"


class SSHProcessAdaptor(ProcessAdaptor, abc.ABC):
    '''The methods of the subinterfaces are all private'''

    ASSETS_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "ssh_process_assets")
    CONVERTERS_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "converters")

    def adapt_job_context(self, job_context):
        process_parameters = job_context["process"]["inputs"]["parameters"]
        print(process_parameters)
        new_process_parameters = {
            SCRIPT_KEY: self.get_script_filenames(),
            SCRIPT_FILES_KEY: self.get_script_files_list(),
            SCRIPT_PARAMS_KEY: self.get_script_params_string(process_parameters),
        }
        job_context["process"]["inputs"]["parameters"] = new_process_parameters
        job_context[RESULTS_KEY] = self.get_results_files_list(process_parameters)
        return job_context

    @abc.abstractmethod
    def get_script_filenames(self):
        raise NotImplementedError

    @abc.abstractmethod
    def get_script_files_list(self):
        raise NotImplementedError

    @abc.abstractmethod
    def get_script_params_string(self, process_parameters):
        raise NotImplementedError

    @abc.abstractmethod
    def get_results_files_list(self, process_parameters):
        raise NotImplementedError
