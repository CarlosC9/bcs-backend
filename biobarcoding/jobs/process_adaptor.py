import abc

from biobarcoding.jobs.ssh_process_adaptors.SSHClustalProcessAdaptor import SSHClustalProcessAdaptor
from biobarcoding.rest import galaxy_tm_processes, ssh_tm_processes
from abc import ABC
import os


class ProcessAdaptor(ABC):
    '''The only public method of Process Adaptors'''

    @abc.abstractmethod
    def adapt_job_context(self, job_context):
        raise NotImplementedError


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
        input_filename = job_context["data"]["input_dataset"]["path"]
        process_parameters = job_context["process_params"]["parameters"]
        new_process_parameters = {
            self.SCRIPT_KEY: self._get_script_filename(),
            self.SCRIPT_FILES_KEY: self._get_script_files_list(),
            self.SCRIPT_PARAMS_KEY: self._get_script_params_string(input_filename, process_parameters),
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

class GalaxyProcessAdaptor(ProcessAdaptor, ABC):
    '''The methods of the subinterfaces are all private'''

    def __init__(self):
        pass


class ProcessAdaptorFactory:

    def get(self, jm_type, process_id):
        process_adaptor = None
        if jm_type == "galaxy":
            process_adaptor = self._get_galaxy_process_adaptor(process_id)
        elif jm_type == "ssh":
            process_adaptor = self._get_ssh_process_adaptor(process_id)

        return process_adaptor

    def _get_ssh_process_adaptor(self, process_id):
        ssh_process_adaptor = None
        if ssh_tm_processes[process_id] == "SSHClustalW":
            ssh_process_adaptor = SSHClustalProcessAdaptor()

        return ssh_process_adaptor

    def _get_galaxy_process_adaptor(self, process_id):
        galaxy_process_adaptor = ""
        if galaxy_tm_processes[process_id] == "klustal-1":
            pass
        # TODO complete with rest of galaxy_tm_processes values

        return galaxy_process_adaptor
