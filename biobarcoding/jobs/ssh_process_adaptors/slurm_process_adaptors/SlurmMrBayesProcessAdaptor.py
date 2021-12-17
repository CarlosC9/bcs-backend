import os

from biobarcoding.jobs.ssh_process_adaptors.slurm_process_adaptors import SlurmProcessAdaptor


class SlurmMrBayesProcessAdaptor(SlurmProcessAdaptor):
    INPUT_FILENAME = "aln.nexus"
    OUTPUT_FILENAME_WITHOUT_EXTENSION = "mrbayes"

    def get_script_filenames(self):
        return [{
            "remote_name": "mb.sh",
            "file": os.path.join(self.ASSETS_FOLDER, "mrbayes_assets", "mb.sh"),
            "type": "sh"
        }]

    def get_script_files_list(self):
        return [
            {
                "remote_name": "mb_batch_template.nex",
                "file": os.path.join(self.ASSETS_FOLDER, "mrbayes_assets", "mb_batch_template.nex"),
                "subprocess": "Mr Bayes",
                "type": "nexus"
            },
            {
                "remote_name": "mb_template_writer.py",
                "file": os.path.join(self.ASSETS_FOLDER, "mrbayes_assets", "mb_template_writer.py"),
                "subprocess": "Mr Bayes",
                "type": "python"
            }
        ]

    def parse_script_params(self, process_parameters):
        mrbayes_parameters = process_parameters["script_parameters"]["Mr Bayes"]
        if mrbayes_parameters['start_phylotree']:
            mrbayes_parameters['input_tree_cmd'] = "\'execute start_phylotree.nexus;\'"
        else:
            mrbayes_parameters['input_tree_cmd'] = "\'[No input tree]\'"
        del mrbayes_parameters['start_phylotree']
        mrbayes_parameters['filename'] = self.OUTPUT_FILENAME_WITHOUT_EXTENSION
        hpc_parameters = process_parameters["hpc_parameters"]
        mrbayes_parameters['ntasks'] = hpc_parameters['ntasks']

        return mrbayes_parameters

    def get_results_files_list(self, process_parameters):
        return [
            {
                "remote_name": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}_phylotrees.nexus",
                "file": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}_phylotrees.nexus",
                "subprocess": "Mr Bayes",
                "object_type": {"bos": "phylotrees"},
                "content_type": "text/x-nexus",
                "type": "nexus"
            },
            {
                "remote_name": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}_substitution_model_parameters.tsv",
                "file": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}_substitution_model_parameters.tsv",
                "subprocess": "Mr Bayes",
                "object_type": {"bos": "mrbayes_substitution_model_parameters"},
                "content_type": "text/tab-separated-values",
                "type": "tsv"
            },
            {
                "remote_name": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}_convergence_diagnostics.tsv",
                "file": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}_convergence_diagnostics.tsv",
                "subprocess": "Mr Bayes",
                "object_type": {"bos": "mrbayes_convergence_diagnostics"},
                "content_type": "text/tab-separated-values",
                "type": "tsv"
            }

        ]
