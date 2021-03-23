from biobarcoding.jobs.galaxy_process_adaptors import  GalaxyProcessAdaptor


job_context = {
  "resource_id": "8fac3ce8-8796-445f-ac27-4baedadeff3b",
  "process_id": "c8df0c20-9cd5-499b-92d4-5fb35b5a369a",
  "process_params": {
    "parameters": {
        {"clustalw":
             {"dnarna": "DNA",
              "outform": "phylip",
              "out_seqnos": "OFF",
              "out_order": "ALIGNED",
              "mode": "complete",
              "seq_range_start":"1",
              "seq_range_end": "99999"}},
    },
    "data":   [
                             {
                         "name": "fasta.fasta",
                         "file" : {
                                {"bo_type": "collection",
                                "ids": ["jkdjdjdjjdjd"],
                                "type": "fasta"},
                                  {
                                    "bo_type": "path",
                                      "ids":["...."],
                                      "type":"fasta"
                                  },
                                 {"bo_type": "seq",
                                  "bo": ["jkdjdjdjjdjd", "jdjdjdjfjdjdj"],
                                "type": "fasta"}
                                }
                             }
                         ]
  },
  "credentials": {
    "api_key": "fakekey"
  }
}

# job_context = {"resource_id":"null",
#                "process_id":"c8df0c20-9cd5-499b-92d4-5fb35b5a369a",
#                "process_params":
#                    {"parameters":
#                         {"clustalw":
#                              {"dnarna":"DNA",
#                               "outform":"phylip",
#                               "out_seqnos":"OFF",
#                               "out_order":"ALIGNED",
#                               "mode":"complete",
#                               "seq_range_start":
#                                   "1","seq_range_end":"99999"}},
#                     "data":
#                         {"Input dataset":
#                              {"path":"dd",
#                               "type":"fasta"}}
#                     },
#                "credentials":
#                    {"api_key":"fakekey"}
#                }

job_context_out = {"endpoint_url": "http//:localhost:5000/",
               "process":
                   {"inputs":
                        {"parameters": "...."},
                         "data": [
                             {
                         "name": "Input dataset",
                         "file" : {
                                {"bo_type": "collection",
                                "ids": ["jkdjdjdjjdjd"],
                                "type": "fasta"},
                                  {
                                    "bo_type": "path",
                                      "ids":["...."],
                                      "type":"fasta"
                                  },
                                 {"bo_type": "seq",
                                  "bo": ["jkdjdjdjjdjd", "jdjdjdjfjdjdj"],
                                "type": "fasta"}
                                }
                             }
                         ]
                    },
                    "name": {"........"},
               "outputs":[{"....."}],
               "status": "created",
               "resource": {"......."},
               "job_id": 60}

class GalaxyClustalAdaptator(GalaxyProcessAdaptor):
    # workflow dict commun to all clustalW
    worflow_dict = {'model_class': 'StoredWorkflow',
 'id': 'df7a1f0c02a5b08e',
 'name': 'MSA ClustalW',
 'create_time': '2021-03-11T14:03:04.319361',
 'update_time': '2021-03-16T08:36:03.759635',
 'published': False,
 'deleted': False,
 'tags': [],
 'latest_workflow_uuid': '7051fb24-072f-4aca-800e-97e5b312a493',
 'url': '/api/workflows/df7a1f0c02a5b08e',
 'owner': 'admin',
 'inputs': {'0': {'label': 'Input dataset',
   'value': '',
   'uuid': 'b2337650-823a-4b85-ab8d-6eaf22601515'}},
 'annotation': 'Multi alignment using ClustalW tool',
 'steps': {'0': {'id': 0,
   'type': 'data_input',
   'tool_id': None,
   'tool_version': None,
   'annotation': None,
   'tool_inputs': {'optional': False},
   'input_steps': {}},
  '1': {'id': 1,
   'type': 'tool',
   'tool_id': 'toolshed.g2.bx.psu.edu/repos/devteam/clustalw/clustalw/2.1',
   'tool_version': '2.1',
   'annotation': None,
   'tool_inputs': {'__page__': None,
    '__rerun_remap_job_id__': None,
    'dnarna': 'DNA',
    'input': {'__class__': 'RuntimeValue'},
    'out_order': 'ALIGNED',
    'outcontrol': {'__current_case__': 2,
     'out_seqnos': 'false',
     'outform': 'clustal'},
    'range': {'__current_case__': 0, 'mode': 'complete'}},
   'input_steps': {'input': {'source_step': 0, 'step_output': 'output'}}}},
 'version': 0}
    def  __complete_inputs_with_labels(self,job_context):
        datas = job_context["process"]["inputs"]["data"]
        data = []
        workspace = str(job_context['job_id'])
        for d in datas:
            data.append(dict(remote_name = d['name'], path = f"{workspace}/", type  = 'fasta', file = d['file']))
        return data

    def __complete_with_outputs_files(self,job_context):# ??
        output_format = job_context['proces']['inputs']['parameters']["clustalw"].get("outform")
        extension = output_format
        output_1 = dict(format = output_format, name = '.'.join([output_format,extension]), path = "" )
        output_2 = dict(format = 'nhx', name = '.'.join(['dnd','nhx']), path = "")
        return [output_1,output_2]

    def __complete_with_errors_files(self, job_context):
        pass
