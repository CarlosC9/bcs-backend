import os

from ..ssh_process_adaptors import SSHProcessAdaptor


class SSHPaupParsimonyProcessAdaptor(SSHProcessAdaptor):
    INPUT_FILENAME = "paup_assets/aln.nexus"

    def get_script_filenames(self):
        return [{
                  "remote_name": "paup_parsimony.sh",
                  "file": os.path.join(self.ASSETS_FOLDER, "paup_parsimony.sh"),
                  "type": "sh"
                }]

    def get_script_files_list(self):
        return [
            {
                "remote_name": "paup_parsimony.nexus",
                "file": os.path.join(self.ASSETS_FOLDER, "paup_assets", "paup_parsimony.nexus"),
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

    def parse_script_params(self, process_parameters):
        parameters = process_parameters["PAUP Parsimony"]
        sets = parameters.get('taxset').sets
        assumptions = "none"
        taxset_paup = parameters.get('taxset').paup
        params_str = f"{parameters.get('outRoot')} {parameters.get('gapMode')} " + \
                     f"{parameters.get('addseq')} {parameters.get('swap')} " + \
                     f"{parameters.get('hold')} \"{parameters.get('consensus_tree_type')}\" " + \
                     f"{parameters.get('le50')} {parameters.get('percent')} " + \
                     f"{parameters.get('search')} {parameters.get('nReplicas')} " + \
                     f"{parameters.get('method')} \"{parameters.get('taxset').enforce_converse}\" " + \
                     f"\"{taxset_paup}\" \"{sets}\" \"{assumptions}\""

        return params_str

    def get_results_files_list(self, process_parameters):
        consensus_tree_filename = ""
        parameters = process_parameters["PAUP Parsimony"]
        files = []
        if parameters.get("method") == "bootstrap":
            consensus_tree_filename = "bootstrap_consensus.nexus"
            files = [
                {
                    "remote_name": f"bootstrap_replicas.nexus",
                    "file": f"bootstrap_replicas.nexus",
                    "subprocess": "PAUP Parsimony",
                    "object_type": {"bos": "phylotrees"},
                    "content_type": "text/x-nexus",
                    "type": "nexus"
                },
            ]
        elif parameters.get("method") == "jackknife":
            consensus_tree_filename = "jackknife_consensus.nexus"
            files = [
                {
                    "remote_name": f"jackknife_replicas.nexus",
                    "file": f"jackknife_replicas.nexus",
                    "subprocess": "PAUP Parsimony",
                    "object_type": {"bos": "phylotrees"},
                    "content_type": "text/x-nexus",
                    "type": "nexus"
                },
            ]
        elif parameters.get("method") == "simple":
            consensus_tree_filename = "hsearch_consensus.nexus"
        return files + [
            {
                "remote_name": consensus_tree_filename,
                "file": consensus_tree_filename,
                "subprocess": "PAUP Parsimony",
                "object_type": {"bos": "phylotrees"},
                "content_type": "text/x-nexus",
                "type": "nexus"
            },
            {
                "remote_name": f"cladogram.txt",
                "file": "cladogram.txt",
                "subprocess": "PAUP Parsimony",
                "object_type": {"bos": "paup_cladogram"},
                "content_type": "text/plain",
                "type": "txt"
            },
            {
                "remote_name": f"treescores.txt",
                "file": "treescores.txt",
                "subprocess": "PAUP Parsimony",
                "object_type": {"bos": "paup_scores"},
                "content_type": "text/tab-separated-values",
                "type": "txt"
            },
            {
                "remote_name": f"ngd_paup_parsimony.nexus",
                "file": "ngd_paup_parsimony.nexus",
                "subprocess": "PAUP Parsimony",
                "object_type": {"bos": "paup_script"},
                "content_type": "text/plain",
                "type": "txt"
            },
            {
                "remote_name": f"sets_and_assumptions.nexus",
                "file": "sets_and_assumptions.nexus",
                "subprocess": "PAUP Parsimony",
                "object_type": {"bos": "phylotrees"},
                "content_type": "text/plain",
                "type": "txt"
            },
        ]
