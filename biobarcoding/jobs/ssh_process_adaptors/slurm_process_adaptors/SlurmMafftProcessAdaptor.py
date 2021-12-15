import os

from biobarcoding.jobs.ssh_process_adaptors.slurm_process_adaptors import SlurmProcessAdaptor


class SlurmMAFFTProcessAdaptor(SlurmProcessAdaptor):
    INPUT_FILENAME = "input.fasta"

    def get_script_filenames(self):
        return [{
            "remote_name": "mafft.sh",
            "file": os.path.join(self.ASSETS_FOLDER, "mafft.sh"),
            "type": "sh"
        }]

    def get_script_files_list(self):
        return []

    def parse_script_params(self, process_parameters):
        mafft_parameters = process_parameters["script_parameters"]["MAFFT"]
        mafft_parameters['input_filename'] = self.INPUT_FILENAME
        mafft_parameters['output_filename'] = self.get_results_files_list(process_parameters)[0]['remote_name']
        hpc_parameters = process_parameters["hpc_parameters"]
        mafft_parameters['cpus_per_task'] = hpc_parameters['cpus_per_task']
        mafft_parameters['ntasks'] = hpc_parameters['ntasks']

        return mafft_parameters

    def get_results_files_list(self, process_parameters):
        return [
            {
                "remote_name": f"mafft.fasta",
                "file": f"mafft.fasta",
                "subprocess": "MAFFT",
                "object_type": {"bos": "alignments"},
                "content_type": "text/x-fasta",
                "type": "fasta"
            }
        ]
