import os

from ..ssh_process_adaptors import SSHProcessAdaptor
from ..ssh_process_adaptors.SSHClustalProcessAdaptor import SSHClustalProcessAdaptor
from ..ssh_process_adaptors.SSHPaupParsimonyProcessAdaptor import SSHPaupParsimonyProcessAdaptor
from ..ssh_process_adaptors.SSHPdaProcessAdaptor import SSHPdaProcessAdaptor


class SSHClustal_PaupParsimony_PdaProcessAdaptor(SSHProcessAdaptor):

    def __init__(self):
        self.clustalAdaptor = SSHClustalProcessAdaptor()
        self.paupParsimonyAdaptor = SSHPaupParsimonyProcessAdaptor()
        self.pdaAdaptor = SSHPdaProcessAdaptor()

    def get_script_filenames(self):
        return [{
            "remote_name": "clustalw+paup_parsimony+pda.sh",
            "file": os.path.join(self.ASSETS_FOLDER, "clustalw+paup_parsimony+pda.sh"),
            "type": "sh"
        }] + self.clustalAdaptor.get_script_filenames() + self.paupParsimonyAdaptor.get_script_filenames() + \
               self.pdaAdaptor.get_script_filenames()

    def get_script_files_list(self):
        return self.clustalAdaptor.get_script_files_list() + \
               self.paupParsimonyAdaptor.get_script_files_list() + \
               self.pdaAdaptor.get_script_files_list() + \
               [{
                   "remote_name": "fasta2nexus.py",
                   "file": os.path.join(self.CONVERTERS_FOLDER, "fasta2nexus.py"),
                   "subprocess": "MSA ClustalW + PAUP Parsimony + Phylogenetic Diversity Analyzer",
                   "type": "python"
               },
                   {
                       "remote_name": "nexus2newick.py",
                       "file": os.path.join(self.CONVERTERS_FOLDER, "nexus2newick.py"),
                       "subprocess": "MSA ClustalW + PAUP Parsimony + Phylogenetic Diversity Analyzer",
                       "type": "python"
                   }]

    def parse_script_params(self, process_parameters):
        return f"{self.clustalAdaptor.parse_script_params(process_parameters)} " + \
               f"{self.paupParsimonyAdaptor.parse_script_params(process_parameters)} " + \
               f"{self.pdaAdaptor.parse_script_params(process_parameters)}"

    def get_results_files_list(self, process_parameters):
        return self.clustalAdaptor.get_results_files_list(process_parameters) + \
               self.paupParsimonyAdaptor.get_results_files_list(process_parameters) + \
               self.pdaAdaptor.get_results_files_list(process_parameters)
