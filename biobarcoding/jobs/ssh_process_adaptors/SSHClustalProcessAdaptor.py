import os

from biobarcoding.jobs.ssh_process_adaptors import SSHProcessAdaptor


class SSHClustalProcessAdaptor(SSHProcessAdaptor):
    INPUT_FILENAME = "input_dataset"

    def get_script_filenames(self):
        return [{
                  "remote_name": "clustalw.sh",
                  "file": os.path.join(self.ASSETS_FOLDER, "clustalw.sh"),
                  "type": "sh"
                }]

    def get_script_files_list(self):
        return []

    def get_script_params_string(self, process_parameters):
        clustalw_parameters = process_parameters["MSA ClustalW"]
        output_file = self.get_results_files_list(process_parameters)[0].get("remote_name")
        params_str = f"{self.INPUT_FILENAME} {output_file} " + \
                     f"{clustalw_parameters['out_order']} {clustalw_parameters['dnarna']} " + \
                     f"{clustalw_parameters['outform']}"
        if clustalw_parameters["mode"] == "part":
            params_str += f" {clustalw_parameters['seq_range_start']} {clustalw_parameters['seq_range_end']}"

        return params_str

    def get_results_files_list(self, process_parameters):
        clustalw_parameters = process_parameters["MSA ClustalW"]
        return [
            {
                "remote_name": f"clustal.{clustalw_parameters['outform'].lower()}",
                "file": f"clustal.{clustalw_parameters['outform'].lower()}",
                "subprocess": "MSA ClustalW",
                "type": "fa"
            },
            {
                "remote_name": f"{self.INPUT_FILENAME}.dnd",
                "file": "dnd.nhx",
                "subprocess": "MSA ClustalW",
                "type": "dnd"
            }
        ]
