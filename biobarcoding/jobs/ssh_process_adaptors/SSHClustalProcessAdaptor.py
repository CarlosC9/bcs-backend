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

    def __complete_input_file_names_list(self, process_params):
        print(self.ASSETS_FOLDER)
        pass

    def __get_script_params_string(self, job_context):
        pass

    def __get_results_file_names_list(self):
        pass




