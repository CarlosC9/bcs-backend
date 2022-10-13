import os

from biobarcoding.jobs.cipres_process_adaptors import CipresProcessAdaptor


class CipresMrBayesProcessAdaptor(CipresProcessAdaptor):
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
                "remote_name": "mb.py",
                "file": os.path.join(self.ASSETS_FOLDER, "mrbayes_assets", "mb.py"),
                "subprocess": "Mr Bayes",
                "type": "python"
            },
            {
                "remote_name": "submit_cipres.py",
                "file": os.path.join(self.ASSETS_FOLDER, "submit_cipres.py"),
                "subprocess": "Mr Bayes",
                "type": "python"
            }
        ]

    def parse_script_params(self, process_parameters):
        mrbayes_parameters = process_parameters["script_params"]["Mr Bayes"]
        '''hpc_parameters = process_parameters["hpc_params"]
        mrbayes_parameters['cpus_per_task'] = hpc_parameters['cpus_per_task']
        mrbayes_parameters['time'] = hpc_parameters['time']'''
        mrbayes_parameters['filename'] = self.OUTPUT_FILENAME_WITHOUT_EXTENSION

        return self.parse_dict_env_variables(mrbayes_parameters)

    def get_results_files_list(self, process_parameters):
        return [
            {
                "remote_name": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}.t",
                "file": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}_phylotrees.nexus",
                "subprocess": "Mr Bayes",
                "object_type": {"bos": "phylotrees"},
                "content_type": "text/x-nexus",
                "type": "nexus"
            },
            {
                "remote_name": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}.p",
                "file": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}_substitution_model_parameters.tsv",
                "subprocess": "Mr Bayes",
                "object_type": {"bos": "mrbayes_substitution_model_parameters"},
                "content_type": "text/tab-separated-values",
                "type": "tsv"
            },
            {
                "remote_name": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}.mcmc",
                "file": f"{self.OUTPUT_FILENAME_WITHOUT_EXTENSION}_convergence_diagnostics.tsv",
                "subprocess": "Mr Bayes",
                "object_type": {"bos": "mrbayes_convergence_diagnostics"},
                "content_type": "text/tab-separated-values",
                "type": "tsv"
            }
        ]
