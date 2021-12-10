import abc
import os

from biobarcoding.jobs.ssh_process_adaptors import SSHProcessAdaptor


class SlurmProcessAdaptor(SSHProcessAdaptor, abc.ABC):

    '''The methods of the subinterfaces are all private'''
    HPC_PARAMETERS_KEY = "hpc_parameters"

    ASSETS_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "slurm_process_assets")
    CONVERTERS_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../ssh_process_assets",
                                     "converters")

    def adapt_job_context(self, job_context):
        script_parameters = job_context["process"]["inputs"]["parameters"]["script_params"]
        hpc_parameters = job_context["process"]["inputs"]["parameters"]["hpc_params"]
        new_process_parameters = {
            self.SCRIPT_KEY: self.get_script_filenames(),
            self.SCRIPT_FILES_KEY: self.get_script_files_list(),
            self.SCRIPT_PARAMS_KEY: {"process_parameters": self.parse_script_params({
                    "script_parameters": script_parameters, "hpc_parameters": hpc_parameters}),
                "hpc_parameters": hpc_parameters},
            self.HPC_PARAMETERS_KEY: hpc_parameters
        }
        job_context["process"]["inputs"]["parameters"] = new_process_parameters
        job_context[self.RESULTS_KEY] = self.get_results_files_list(script_parameters)
        return job_context
