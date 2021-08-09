import os

from ..ssh_process_adaptors import SSHProcessAdaptor
from ..ssh_process_adaptors.SSHPaupParsimonyProcessAdaptor import SSHPaupParsimonyProcessAdaptor
from ..ssh_process_adaptors.SSHPdaProcessAdaptor import SSHPdaProcessAdaptor


class SSHPaupParsimony_PdaProcessAdaptor(SSHProcessAdaptor):

    def __init__(self):
        self.paupParsimonyAdaptor = SSHPaupParsimonyProcessAdaptor()
        self.pdaAdaptor = SSHPdaProcessAdaptor()

    def get_script_filenames(self):
        return [{
            "remote_name": "paup_parsimony+pda.sh",
            "file": os.path.join(self.ASSETS_FOLDER, "paup_parsimony+pda.sh"),
            "type": "sh"
        }] + self.paupParsimonyAdaptor.get_script_filenames() + self.pdaAdaptor.get_script_filenames()

    def get_script_files_list(self):
        return self.paupParsimonyAdaptor.get_script_files_list() + \
               self.pdaAdaptor.get_script_files_list() + \
               [{
                   "remote_name": "nexus2newick.py",
                   "file": os.path.join(self.CONVERTERS_FOLDER, "nexus2newick.py"),
                   "subprocess": "PAUP Parsimony + Phylogenetic Diversity Analyzer",
                   "type": "txt"
               }]

    def get_script_params_string(self, process_parameters):
        return f"{self.paupParsimonyAdaptor.get_script_params_string(process_parameters)} " + \
               f"{self.pdaAdaptor.get_script_params_string(process_parameters)}"

    def get_results_files_list(self, process_parameters):
        return self.paupParsimonyAdaptor.get_results_files_list(process_parameters) + \
               self.pdaAdaptor.get_results_files_list(process_parameters)
