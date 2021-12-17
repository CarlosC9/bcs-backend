import os

from biobarcoding.jobs.ssh_process_adaptors import SSHProcessAdaptor


class SSHPdaProcessAdaptor(SSHProcessAdaptor):

    def get_script_filenames(self):
        return [{
            "remote_name": "pda.sh",
            "file": os.path.join(self.ASSETS_FOLDER, "pda.sh"),
            "type": "sh"
        }]

    def get_script_files_list(self):
        return [
            {
                "remote_name": "consensus.newick",
                "file": os.path.join(self.ASSETS_FOLDER, "pda_assets", "heuristic_phylogeny.newick"),
                "subprocess": "Phylogenetic Diversity Analyzer",
                "type": "newick"
            }
        ]

    def parse_script_params(self, process_parameters):
        pda_parameters = process_parameters["Phylogenetic Diversity Analyzer"]
        params_str = ""
        params_str += "-endem " if pda_parameters['endem'] else params_str
        params_str += "-excl" if pda_parameters['excl'] else params_str
        return params_str

    def get_results_files_list(self, process_parameters):
        return [
            {
                "remote_name": "output.pda",
                "file": "output.pda",
                "subprocess": "Phylogenetic Diversity Analyzer",
                "object_type": {"geo": "layers"},
                "content_type": "text/plain+geolayer",
                "type": "pda"
            }
        ]
