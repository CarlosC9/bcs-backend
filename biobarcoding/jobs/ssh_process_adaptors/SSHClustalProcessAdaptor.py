from biobarcoding.jobs.ssh_process_adaptors import SSHProcessAdaptor


class SSHClustalProcessAdaptor(SSHProcessAdaptor):
    INPUT_FILENAME = "input_dataset"

    def _get_script_filename(self):
        return "clustalw.sh"

    def _get_script_files_list(self):
        return []

    def _get_script_params_string(self, process_parameters):
        clustalw_parameters = process_parameters["MSA ClustalW"]
        output_file = self._get_results_files_list(process_parameters)[0].get("remote_name")
        params_str = f"{self.INPUT_FILENAME} {output_file} " + \
                     f"{clustalw_parameters['out_order']} {clustalw_parameters['dnarna']} " + \
                     f"{clustalw_parameters['outform']}"
        if clustalw_parameters["mode"] == "part":
            params_str += f" {clustalw_parameters['seq_range_start']} {clustalw_parameters['seq_range_end']}"

        return params_str

    def _get_results_files_list(self, process_parameters):
        clustalw_parameters = process_parameters["MSA ClustalW"]
        return [
            {
                "remote_name": f"clustal.{clustalw_parameters['outform'].lower()}",
                "file": f"clustal.{clustalw_parameters['outform'].lower()}",
                "type": "fa"
            },
            {
                "remote_name": f"{self.INPUT_FILENAME}.dnd",
                "file": "dnd.nhx",
                "type": "dnd"
            }
        ]
