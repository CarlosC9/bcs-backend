import os
from biobarcoding.jobs.ssh_process_adaptors import SSHProcessAdaptor

class SSHClustalProcessAdaptor(SSHProcessAdaptor):

    INPUT_FILENAME = "clustalw.fasta"

    def _get_script_filename(self):
        return "clustalw.sh"

    def _get_script_files_list(self):
        return [
            {
                "remote_name": "clustalw.sh",
                "file": os.path.join(self.ASSETS_FOLDER, "clustalw.sh"),
                "type": "sh"
            }
        ]

    def _get_script_params_string(self, process_parameters):
        output_file = self._get_results_files_list()[0].get("remote_name")
        params_str = f"{self.INPUT_FILENAME} {output_file} " +\
                     f"{process_parameters['out_order']} {process_parameters['dnarna']}"
        if process_parameters["mode"] == "part":
            params_str += f" {process_parameters['seq_range_start']} {process_parameters['seq_range_end']}"

        return params_str

    def _get_results_files_list(self):
        return [
            {
                "remote_name": "clustalw.aln",
                "file": "clustalw.aln",
                "type": "aln"
            },
            {
                "remote_name": "clustalw.dnd",
                "file": "clustalw.dnd",
                "type": "dnd"
            }
        ]

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
                        "remote_name": "clustalw.aln",
                        "file": "clustalw.aln",
                        "type": "aln"
                    },
                    {
                        "remote_name": "clustalw.dnd",
                        "file": "clustalw.dnd",
                        "type": "dnd"
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



