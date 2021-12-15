import os

from ..ssh_process_adaptors import SSHProcessAdaptor


class SSHClustalProcessAdaptor(SSHProcessAdaptor):
    INPUT_FILENAME = "input_dataset.fasta"

    def get_script_filenames(self):
        return [{
            "remote_name": "clustalw.sh",
            "file": os.path.join(self.ASSETS_FOLDER, "clustalw.sh"),
            "type": "sh"
        }]

    def get_script_files_list(self):
        return []

    def parse_script_params(self, process_parameters):
        clustalw_parameters = process_parameters["MSA ClustalW"]
        output_file = self.get_results_files_list(process_parameters)[0].get("remote_name")
        params_str = f"{self.INPUT_FILENAME} {output_file} " + \
                     f"{clustalw_parameters['out_order']} {clustalw_parameters['dnarna']} " + \
                     f"{clustalw_parameters['outform']}"
        if clustalw_parameters["mode"] == "part":
            params_str += f" {clustalw_parameters['seq_range_start']} {clustalw_parameters['seq_range_end']}"
        else:
            params_str += " none none"

        return params_str

    def get_results_files_list(self, process_parameters):
        clustalw_parameters = process_parameters["MSA ClustalW"]
        return [
            {
                "remote_name": f"clustal.{clustalw_parameters['outform'].lower()}",
                "file": f"clustal.{clustalw_parameters['outform'].lower()}",
                "subprocess": "MSA ClustalW",
                "object_type": {"bos": "alignments"},
                "content_type": f"text/x-{clustalw_parameters['outform'].lower()}",
                "type": clustalw_parameters['outform'].lower()
            },
            {
                "remote_name": f"input_dataset.newick",
                "file": "intermidiate_tree.newick",
                "subprocess": "MSA ClustalW",
                "object_type": {"bos": "phylotrees"},
                "content_type": "text/x-nhx",
                "type": "newick"
            }
        ]
