import abc

from biobarcoding.jobs.ssh_process_adaptors.SSHTestProcessAdaptor import SSHTestProcessAdaptor
from biobarcoding.rest import galaxy_tm_processes, ssh_tm_processes
from abc import ABC


class ProcessAdaptor(ABC):
    '''The only public method of Process Adaptors'''

    @abc.abstractmethod
    def fill_job_context_with_process(self, job_context):
        raise NotImplementedError


class SSHProcessAdaptor(ProcessAdaptor, ABC):
    '''The methods of the subinterfaces are all private'''

    @abc.abstractmethod
    def __complete_input_file_names_list(self, job_context):
        raise NotImplementedError

    @abc.abstractmethod
    def __get_script_params_string(self, job_context):
        raise NotImplementedError

    @abc.abstractmethod
    def __get_results_file_names_list(self):
        raise NotImplementedError


class GalaxyProcessAdaptor(ProcessAdaptor, ABC):
    '''The methods of the subinterfaces are all private'''

    def __init__(self):
        pass


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
            ssh_process_adaptor = SSHTestProcessAdaptor()

        return ssh_process_adaptor

    def get_galaxy_process_adaptor(self, process_id):
        galaxy_process_adaptor = None
        if galaxy_tm_processes[process_id] == "klustal-1":
            pass
        # TODO complete with rest of galaxy_tm_processes values

        return galaxy_process_adaptor
