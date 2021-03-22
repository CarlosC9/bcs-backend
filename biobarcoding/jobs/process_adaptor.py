import abc

from biobarcoding.jobs.ssh_process_adaptors.SSHClustalProcessAdaptor import SSHClustalProcessAdaptor
from biobarcoding.jobs.galaxy_process_adaptors.GalaxyClustalProcessAdaptor import GalaxyClustalAdaptator
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
        process_parameters = job_context["process_params"]["parameters"]
        new_process_parameters = {
            self.SCRIPT_KEY: self.__get_script_filename(),
            self.SCRIPT_FILES_KEY: self.__get_script_files_list(),
            self.SCRIPT_PARAMS_KEY: self.__get_script_params_string(process_parameters),
            self.RESULTS_FILES_KEY: self.__get_results_files_list(),
        }
        job_context["process_params"]["parameters"] = new_process_parameters
        return job_context

    @abc.abstractmethod
    def __get_script_filename(self):
        raise NotImplementedError

    @abc.abstractmethod
    def __get_script_files_list(self):
        raise NotImplementedError

    @abc.abstractmethod
    def __get_script_params_string(self, process_parameters):
        raise NotImplementedError

    @abc.abstractmethod
    def __get_results_files_list(self):
        raise NotImplementedError

class GalaxyProcessAdaptor(ProcessAdaptor, ABC):
    '''The methods of the subinterfaces are all private'''

    @abc.abstractmethod
    def __complete_inputs_with_labels(self,job_context):
        raise NotImplementedError
    @abc.abstractmethod
    def __complete_with_outputs_files(self,job_context):
        raise NotImplementedError

    def __complete_with_errors_files(self, job_context):
        raise NotImplementedError

    def adapt_job_context(self, job_context):
        job_context['process']['inputs']['data'] = self.__complete_inputs_with_labels(job_context)
        job_context['results'] = self.__complete_with_outputs_files(job_context)
        return job_context


class ProcessAdaptorFactory:

    def get(self, jm_type, process_id):
        process_adaptor = None
        if jm_type == "galaxy":
            process_adaptor = self.getGalaxyProcessAdaptor(process_id)
        elif jm_type == "ssh":
            process_adaptor = self.getSSHProcessAdaptor(process_id)

        return process_adaptor

    def get_ssh_process_adaptor(self, process_id):
        ssh_process_adaptor = None
        if ssh_tm_processes[process_id] == "SSHTestProcess":
            ssh_process_adaptor = SSHClustalProcessAdaptor()

        return ssh_process_adaptor

    def get_galaxy_process_adaptor(self, process_id):
        galaxy_process_adaptor = None
        if galaxy_tm_processes[process_id] == "MSA ClustalW":
            galaxy_process_adaptor = GalaxyClustalAdaptator()
        return galaxy_process_adaptor
