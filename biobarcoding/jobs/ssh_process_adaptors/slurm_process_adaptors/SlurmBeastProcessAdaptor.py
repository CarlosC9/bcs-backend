import os

from biobarcoding.jobs.ssh_process_adaptors.slurm_process_adaptors import SlurmProcessAdaptor


class SlurmBeastProcessAdaptor(SlurmProcessAdaptor):
    OUTPUT_FILENAME_WITHOUT_EXTENSION = "mrbayes"

    def get_script_filenames(self):
        return [{
            "remote_name": "beast.sh",
            "file": os.path.join(self.ASSETS_FOLDER, "beast_assets", "beast.sh"),
            "type": "sh"
        }]

    def get_script_files_list(self):
        return [
            {
                "remote_name": "beast_loop.sh",
                "file": os.path.join(self.ASSETS_FOLDER, "beast_assets", "beast_loop.sh"),
                "subprocess": "Beast",
                "type": "sh"
            },
            {
                "remote_name": "beast_template.bea",
                "file": os.path.join(self.ASSETS_FOLDER, "beast_assets", "beast_template.bea"),
                "subprocess": "Beast",
                "type": "beauti"
            },
            {
                "remote_name": "complete_beasy_template.py",
                "file": os.path.join(self.ASSETS_FOLDER, "beast_assets", "complete_beasy_template.py"),
                "subprocess": "Beast",
                "type": "python"
            },
            {
                "remote_name": "create_bayes_factor.py",
                "file": os.path.join(self.ASSETS_FOLDER, "beast_assets", "create_bayes_factor.py"),
                "subprocess": "Beast",
                "type": "python"
            },
            {
                "remote_name": "create_path_sampling.py",
                "file": os.path.join(self.ASSETS_FOLDER, "beast_assets", "create_path_sampling.py"),
                "subprocess": "Beast",
                "type": "python"
            },
            {
                "remote_name": "estimateClockRates.py",
                "file": os.path.join(self.ASSETS_FOLDER, "beast_assets", "estimateClockRates.py"),
                "subprocess": "Beast",
                "type": "python"
            },
            {
                "remote_name": "fasta2nexus.py",
                "file": os.path.join(self.CONVERTERS_FOLDER, "fasta2nexus.py"),
                "subprocess": "Beast",
                "type": "python"
            },
            {
                "remote_name": "modeltest.py",
                "file": os.path.join(self.ASSETS_FOLDER, "beast_assets", "modeltest.py"),
                "subprocess": "Beast",
                "type": "python"
            },
            {
                "remote_name": "nexus2newick_beast_annotations.py",
                "file": os.path.join(self.ASSETS_FOLDER, "beast_assets", "nexus2newick_beast_annotations.py"),
                "subprocess": "Beast",
                "type": "python"
            },
            {
                "remote_name": "tracerer.R",
                "file": os.path.join(self.ASSETS_FOLDER, "beast_assets", "tracerer.R"),
                "subprocess": "Beast",
                "type": "R"
            },
        ]

    def parse_script_params(self, process_parameters):
        beast_parameters = process_parameters["script_parameters"]["Beast"]
        beast_parameters['alignments'] = {}
        parsed_monophyly = {}
        for key, value in beast_parameters['monophyly'].items():
            parsed_monophyly[key] = " ".join(value)
        beast_parameters['monophyly'] = parsed_monophyly
        for d in process_parameters['data']:
            if 'region' in d.keys():
                beast_parameters['alignments'][d['remote_name'].split('.')[0]] = d['region']
        hpc_parameters = process_parameters["hpc_parameters"]
        beast_parameters['threads'] = hpc_parameters['cpus_per_task']

        beast_env_vars = self.parse_dict_env_variables(beast_parameters)
        return beast_env_vars

    def get_results_files_list(self, process_parameters):
        return [
            {
                "remote_name": f"coupled_mcmc/BirthDeathModel/birthDeathModel.nexus",
                "file": f"BirthDeathModel.nwk",
                "subprocess": "Beast",
                "object_type": {"bos": "phylotrees"},
                "content_type": "text/newick",
                "type": "nexus"
            },
            {
                "remote_name": f"coupled_mcmc/YuleModel/yuleModel.nexus",
                "file": f"YuleModel.nwk",
                "subprocess": "Beast",
                "object_type": {"bos": "phylotrees"},
                "content_type": "text/tab-separated-values",
                "type": "nexus"
            },
            {
                "remote_name": f"bayes_factor.txt",
                "file": f"bayes_factor.txt",
                "subprocess": "Beast",
                "object_type": {"bos": "bayes_factor"},
                "content_type": "text",
                "type": "txt"
            }

        ]
