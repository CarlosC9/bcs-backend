import abc
import os

from ..process_adaptor import ProcessAdaptor


class CipresProcessAdaptor(ProcessAdaptor, abc.ABC):
    SCRIPT_KEY = "scripts"
    SCRIPT_FILES_KEY = "scripts_files"
    SCRIPT_PARAMS_KEY = "script_params"
    RESULTS_KEY = "results"

    '''The methods of the subinterfaces are all private'''

    ASSETS_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "cipres_process_assets")

    def adapt_job_context(self, job_context):
        process_parameters = job_context["process"]["inputs"]["parameters"]["script_params"]
        print(process_parameters)
        new_process_parameters = {
            self.SCRIPT_KEY: self.get_script_filenames(),
            self.SCRIPT_FILES_KEY: self.get_script_files_list(),
            self.SCRIPT_PARAMS_KEY: self.parse_script_params(process_parameters),
        }
        job_context["process"]["inputs"]["adapted_parameters"] = new_process_parameters
        job_context[self.RESULTS_KEY] = self.get_results_files_list(process_parameters)
        return job_context

    @abc.abstractmethod
    def get_script_filenames(self):
        raise NotImplementedError

    @abc.abstractmethod
    def get_script_files_list(self):
        raise NotImplementedError

    @abc.abstractmethod
    def parse_script_params(self, process_parameters):
        raise NotImplementedError

    @abc.abstractmethod
    def get_results_files_list(self, process_parameters):
        raise NotImplementedError
