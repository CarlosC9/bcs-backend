import os
from biobarcoding.jobs.process_adaptor import SSHProcessAdaptor

class SSHClustalProcessAdaptor(SSHProcessAdaptor):
    """
    {
       "resource_id":"0292821a-dd33-450a-bdd8-813b2b95c456",
       "process_id":"25932546-d26c-4367-8c81-0c682094d117",
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
          "data":[
             {
                "remote_name": clustalw.fasta,
                "file": [
                    {
                        "bo_type": "collection",
                        "ids": ["jkdjdjdjjdjd"]
                    },
                    {
                        "bo_type": "path",
                        "ids":["...."]
                    },
                    {
                        "bo_type": "seq",
                        "bo": ["jkdjdjdjjdjd", "jdjdjdjfjdjdj"]
                    }
                ],
                "type":"fasta"
             }
          ]
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
            "remote_name": "clustalw.aln",
            "type": "aln"
        }]




