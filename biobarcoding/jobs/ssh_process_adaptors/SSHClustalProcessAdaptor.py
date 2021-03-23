import os
from biobarcoding.jobs.process_adaptor import SSHProcessAdaptor

class SSHClustalProcessAdaptor(SSHProcessAdaptor):
    """
    {
       "resource_id":"0292821a-dd33-450a-bdd8-813b2b95c456",
       "process_id":"c8df0c20-9cd5-499b-92d4-5fb35b5a369a",
       "process_params":{
          "parameters":{
             "dnarna":"DNA",
             "outform":"fasta",
             "out_seqnos":"false",
             "out_order":"ALIGNED",
             "mode":"complete",
             "seq_range_start":"1",
             "seq_range_end":"99999"
          },
          "data":{
             "input_dataset":{
                "type":"fasta"
             }
          }
       },
       "credentials":{
          "known_hosts_filepath":"/home/daniel/.ssh/known_hosts",
          "username":"dreyes"
       }
    }
    """

    INPUT_FILENAME = "clustalw.fasta"

    def _get_script_filename(self):
        return os.path.join(self.ASSETS_FOLDER, "clustal.sh")

    def _get_script_files_list(self):
        return []

    def _get_script_params_string(self, process_parameters):
        output_file = self._get_results_files_list()[0].get("name")
        params_str = f"${self.INPUT_FILENAME} ${output_file} " +\
                     f"${process_parameters['out_order']} ${process_parameters['dnarna']}"
        if process_parameters["mode"] == "part":
            params_str += f" ${process_parameters['seq_range_start']} ${process_parameters['seq_range_end']}"

    def _get_results_files_list(self, data_list):
        [{
            "name": "clustalw.aln",
            "type": "aln"
        }]




