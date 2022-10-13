import os

from biobarcoding.jobs.cipres_process_adaptors import CipresProcessAdaptor


class CipresMafftProcessAdaptor(CipresProcessAdaptor):
    INPUT_FILENAME = "input.fasta"

    def get_script_filenames(self):
        return [{
            "remote_name": "mafft.sh",
            "file": os.path.join(self.ASSETS_FOLDER, "mafft_assets", "mafft.sh"),
            "type": "sh"
        }]

    def get_script_files_list(self):
        return [
            {
                "remote_name": "mafft.py",
                "file": os.path.join(self.ASSETS_FOLDER, "mafft_assets", "mafft.py"),
                "subprocess": "MAFFT",
                "type": "python"
            },
            {
                "remote_name": "submit_cipres.py",
                "file": os.path.join(self.ASSETS_FOLDER, "submit_cipres.py"),
                "subprocess": "MAFFT",
                "type": "python"
            }
        ]

    def parse_script_params(self, process_parameters):
        mafft_parameters = process_parameters["script_params"]["MAFFT"]
        mafft_parameters['input_filename'] = self.INPUT_FILENAME
        '''hpc_parameters = process_parameters["hpc_params"]
        mafft_parameters['cpus_per_task'] = hpc_parameters['cpus_per_task']
        mafft_parameters['ntasks'] = hpc_parameters['ntasks']
        mafft_parameters['time'] = hpc_parameters['time']'''

        return self.parse_dict_env_variables(mafft_parameters)

    def get_results_files_list(self, process_parameters):
        return [
            {
                "remote_name": f"output.mafft",
                "file": f"mafft.fasta",
                "subprocess": "MAFFT",
                "object_type": {"bos": "alignments"},
                "content_type": "text/x-fasta",
                "type": "fasta"
            }
        ]
