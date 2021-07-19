from biobarcoding.jobs.galaxy_process_adaptors import GalaxyProcessAdaptor

# clustal process job context example

class GalaxyClustalAdaptator(GalaxyProcessAdaptor):
    # workflow dict commun to all clustalW
    INPUT_LABEL = 'input_dataset'
    INPUT_EXTENSION = 'fasta'

    def _complete_inputs_with_labels(self, job_context):
        for data in job_context['process']['inputs']['data']:
            data['remote_name'] = self.INPUT_LABEL
            data['file'] = '.'.join([self.INPUT_LABEL, self.INPUT_EXTENSION])
        return job_context

    def _complete_with_outputs_files(self, job_context):
        extension = job_context['process']['inputs']['parameters']["MSA ClustalW"].get("outform")
        output_1 = dict(type='fa', remote_name=f"ClustalW on data 1: clustal",
                        file='.'.join(["clustal", extension]), subrocess="MSA ClustalW")
        output_2 = dict(type='nhx', remote_name=f"ClustalW on data 1: dnd",
                        file='.'.join(['dnd', 'nhx']), subrocess="MSA ClustalW")  # TODO generalizar con regex
        return [output_1, output_2]

    def _complete_with_worlflow_name(self):
        name = 'MSA ClustalW'
        return name

'''
job_context = {'process':
                   {'inputs':
                        {'parameters':
                             {'dnarna': 'DNA',
                              'outform': 'fasta',
                              'out_seqnos': 'false',
                              'out_order': 'ALIGNED',
                              'mode': 'complete',
                              'seq_range_start': '1',
                              'seq_range_end': '99999'},
                         'data': [
                             {'type': 'fasta',
                              'selection': [1, 2],
                              'bo_type': 'sequences'}
                         ]
                         },
                    'name': 'MSA ClustalW'
                    },
               'status': 'created',
               'endpoint_url': 'http://localhost:5000',
               'resource': {
                   'name': 'localhost - galaxy',
                   'jm_type': 'galaxy',
                   'jm_location':
                       {'url': 'http://localhost:8080/'},
                   'jm_credentials': {'api_key': 'fakekey'}
               },
               'job_id': 3}
'''