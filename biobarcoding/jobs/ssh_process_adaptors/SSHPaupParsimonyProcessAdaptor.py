import os

from dotted.collection import DottedDict

from biobarcoding.jobs.ssh_process_adaptors import SSHProcessAdaptor
import json


class SSHPaupParsimonyProcessAdaptor(SSHProcessAdaptor):
    INPUT_FILENAME = "paup_assets/alignment.txt"

    def _get_script_filename(self):
        return "paup_parsimony.sh"

    def _get_script_files_list(self):
        return [
            {
                "remote_name": "paup_parsimony.txt",
                "file": os.path.join(self.ASSETS_FOLDER, "paup_assets", "paup_parsimony.txt"),
                "subprocess": "PAUP Parsimony",
                "type": "txt"
            },
            {
                "remote_name": "alignment.txt",
                "file": os.path.join(self.ASSETS_FOLDER, "paup_assets", "alignment.txt"),
                "subprocess": "PAUP Parsimony",
                "type": "txt"
            },
            {
                "remote_name": "paup_parsimony_params.py",
                "file": os.path.join(self.ASSETS_FOLDER, "paup_assets", "paup_parsimony_params.py"),
                "subprocess": "PAUP Parsimony",
                "type": "py"
            }
        ]

    def _get_script_params_string(self, process_parameters):
        parameters = process_parameters["PAUP Parsimony"]
        taxset = {'sets': '\t\ttaxset myOutgroup = DROMTTGNC Dros;\n', 'assumptions': ''}
        sets = parameters.get('charset').sets + taxset['sets']
        assumptions = parameters.get('charset').assumptions + taxset['assumptions']
        params_str = f"{parameters.get('outRoot')} {parameters.get('gapMode')} " + \
                     f"{parameters.get('addseq')} {parameters.get('swap')} " + \
                     f"{parameters.get('hold')} \"{parameters.get('consensus_tree_type')}\" " + \
                     f"{parameters.get('le50')} {parameters.get('percent')} " + \
                     f"{parameters.get('search')} {parameters.get('nReplicas')} " + \
                     f"{parameters.get('method')} \"{sets}\" \"{assumptions}\""

        return params_str

    def _get_results_files_list(self, process_parameters):
        consensus_tree_filename = ""
        parameters = process_parameters["PAUP Parsimony"]
        files = []
        if parameters.get("method") == "bootstrap":
            consensus_tree_filename = "bootstrap_consensus.tre"
            files = [
                {
                    "remote_name": f"bootstrap_replicas.tre",
                    "file": f"bootstrap_replicas.tre",
                    "subprocess": "PAUP Parsimony",
                    "type": "nexus"
                },
            ]
        elif parameters.get("method") == "jackknife":
            consensus_tree_filename = "jackknife_consensus.tre"
            files = [
                {
                    "remote_name": f"jackknife_replicas.tre",
                    "file": f"jackknife_replicas.tre",
                    "subprocess": "PAUP Parsimony",
                    "type": "nexus"
                },
            ]
        elif parameters.get("method") == "simple":
            consensus_tree_filename = "hsearch_consensus.tre"
        return files + [
            {
                "remote_name": consensus_tree_filename,
                "file": consensus_tree_filename,
                "subprocess": "PAUP Parsimony",
                "type": "nexus"
            },
            {
                "remote_name": f"treedescription.txt",
                "file": "treedescription.txt",
                "subprocess": "PAUP Parsimony",
                "type": "txt"
            },
            {
                "remote_name": f"treescores.txt",
                "file": "treescores.txt",
                "subprocess": "PAUP Parsimony",
                "type": "txt"
            },
            {
                "remote_name": f"ngd_paup_parsimony.txt",
                "file": "ngd_paup_parsimony.txt",
                "subprocess": "PAUP Parsimony",
                "type": "txt"
            },
            {
                "remote_name": f"sets_and_assumptions.txt",
                "file": "sets_and_assumptions.txt",
                "subprocess": "PAUP Parsimony",
                "type": "txt"
            },
        ]
