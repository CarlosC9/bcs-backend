import os

from biobarcoding.jobs.ssh_process_adaptors.slurm_process_adaptors import SlurmProcessAdaptor


class SlurmBlastProcessAdaptor(SlurmProcessAdaptor):
    INPUT_FILENAME = "query.fasta"

    def get_script_filenames(self):
        return [{
            "remote_name": "blast.sh",
            "file": os.path.join(self.ASSETS_FOLDER, "blast.sh"),
            "type": "sh"
        }]

    def get_script_files_list(self):
        return []

    def parse_script_params(self, process_parameters):
        blast_parameters = process_parameters["script_parameters"]["BLAST"]
        blast_parameters['query_filename'] = self.INPUT_FILENAME
        blast_parameters['output_filename'] = self.get_results_files_list(process_parameters)[0]['remote_name']
        hpc_parameters = process_parameters["hpc_parameters"]
        blast_parameters['cpus_per_task'] = hpc_parameters['cpus_per_task']

        return blast_parameters

    def get_results_files_list(self, process_parameters):
        return [
            {
                "remote_name": f"blast.txt",
                "file": f"blast.txt",
                "subprocess": "BLAST",
                "object_type": {"bos": "blast"},
                "content_type": "text/plain",
                "type": "txt"
            }
        ]
