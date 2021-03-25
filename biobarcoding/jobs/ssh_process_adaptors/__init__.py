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
        new_process_parameters = {
            SCRIPT_KEY: self._get_script_filename(),
            SCRIPT_FILES_KEY: self._get_script_files_list(),
            SCRIPT_PARAMS_KEY: self._get_script_params_string(process_parameters),
            RESULTS_FILES_KEY: self._get_results_files_list(),
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
    def _get_results_files_list(self):
        raise NotImplementedError


"""
{
    "process": {
        "inputs": {
            "parameters": {
                "script": "/home/daniel/Documentos/GIT/bcs-backend/biobarcoding/jobs/ssh_process_adaptors/ssh_process_assets/clustal.sh",
                "script_files": [],
                "script_params": "$clustalw.fasta $None $ALIGNED $DNA",
                "result_files": [
                    {
                        "remote_path": "clustalw.aln",
                        "type": "aln"
                    }
                ]
            }, 
            "data": [
                {
                    "remote_name": "clustalw.fasta",
                    "file": "/home/daniel/Documentos/GIT/bcs-backend/tests/ssh_data_test/clustalw.fasta",
                    "type": "fasta"
                }
            ]
        },
        "name": "MSA ClustalW"
    }, 
    "status": "created",
    "endpoint_url": "http://localhost:5000",
    "resource": {
        "name": "balder - ssh",
        "jm_type": "ssh",
        "jm_location": {"host": "balder"},
        "jm_credentials": {
            "known_hosts_filepath": "/home/daniel/.ssh/known_hosts",
            "username": "dreyes"
        }
    },
    "job_id": 12
}"""
