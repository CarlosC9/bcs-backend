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
        output_file = self._get_results_files_list()[0].get("remote_path")
        params_str = f"{self.INPUT_FILENAME} {output_file} " +\
                     f"{process_parameters['out_order']} {process_parameters['dnarna']}"
        if process_parameters["mode"] == "part":
            params_str += f" {process_parameters['seq_range_start']} {process_parameters['seq_range_end']}"

        return params_str

    def _get_results_files_list(self):
        return [{
            "remote_path": "clustalw.aln",
            "type": "aln"
        }]




