import os

from ..ssh_process_adaptors import SSHProcessAdaptor
from ..ssh_process_adaptors.SSHClustalProcessAdaptor import SSHClustalProcessAdaptor
from ..ssh_process_adaptors.SSHPaupParsimonyProcessAdaptor import SSHPaupParsimonyProcessAdaptor


class SSHClustal_PaupParsimonyProcessAdaptor(SSHProcessAdaptor):
    INPUT_FILENAME = "paup_assets/alignment.txt"

    def __init__(self):
        self.clustalAdaptor = SSHClustalProcessAdaptor()
        self.paupParsimonyAdaptor = SSHPaupParsimonyProcessAdaptor()

    def get_script_filenames(self):
        return [{
            "remote_name": "clustalw+paup_parsimony.sh",
            "file": os.path.join(self.ASSETS_FOLDER, "clustalw+paup_parsimony.sh"),
            "type": "sh"
        }] + self.clustalAdaptor.get_script_filenames() + self.paupParsimonyAdaptor.get_script_filenames()

    def get_script_files_list(self):
        return self.clustalAdaptor.get_script_files_list() + \
               self.paupParsimonyAdaptor.get_script_files_list() + \
               [{
                   "remote_name": "fasta2nexus.py",
                   "file": os.path.join(self.CONVERTERS_FOLDER, "fasta2nexus.py"),
                   "subprocess": "MSA ClustalW + PAUP Parsimony",
                   "type": "txt"
               }]

    def get_script_params_string(self, process_parameters):
        return f"{self.clustalAdaptor.get_script_params_string(process_parameters)} " + \
               f"{self.paupParsimonyAdaptor.get_script_params_string(process_parameters)}"

    def get_results_files_list(self, process_parameters):
        return self.clustalAdaptor.get_results_files_list(process_parameters) + \
               self.paupParsimonyAdaptor.get_results_files_list(process_parameters)
