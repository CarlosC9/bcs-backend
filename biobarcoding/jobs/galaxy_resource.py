from bioblend import galaxy
import time
import yaml
import os
from bioblend.galaxy import workflows

def galaxy_instance(path, name ='__default'):
    data = read_yaml_file(path)
    assert name in data, 'unknown instance'
    gal = data[name]
    if isinstance(gal,dict):
        return gal
    else:
        return data[gal] #caso de que es el defaoult




def login(api_key:'str', url:'str'):
    user_key = api_key
    gi = galaxy.GalaxyInstance(url=url, key=user_key)
    return gi


def library_list(gi):
    libraries = gi.libraries.get_libraries()
    return libraries


def library_id(gi,libname: 'str'):
    libraries = gi.libraries.get_libraries()
    lib = next(item for item in libraries if item["name"] == libname)
    return lib['id']


def workflow_list(gi):
    workflows = gi.workflows.get_workflows()
    return workflows


def workflow_id(gi, name: 'str'):
    workflows = gi.workflows.get_workflows()
    work = next(item for item in workflows if item["name"] == name)
    return work['id']


def get_workflow_from_file(gi, workflow_file):
    import_workflow = [gi.workflows.import_workflow_from_local_path(file_local_path=workflow_file)]
    return import_workflow


def get_workflow_from_name(gi, workflow_name):
    wf = gi.workflows.get_workflows(name=workflow_name)
    return wf


def get_workflow_id(wf):
    for wf_dic in wf:
        wf_id = wf_dic['id']
    return wf_id

def read_yaml_file(yaml_path):
    """
    Reads a YAML file safely.

    :param yaml_path:
    :return: dictionary object from YAML content.
    """
    with open(yaml_path, 'r') as f:
        stream = f.read()
    return yaml.safe_load(stream)

def dataset_list(gi):
    datasets = gi.datasets.get_datasets()
    return datasets


def dataset_id(gi, name: 'str'):
    datasets = dataset_list(gi)
    dat = next(item for item in datasets if item["name"] == name)
    return dat['id']

def parse_input(input: 'str', ext : 'str'):
    '''
    :param gi:
    :param input:
    :return:
    '''
    str = input.split(".")
    if str[1] != ext:
        return 'wrong input file'

def create_library(gi, library_name: 'str'):
    '''
    look for a library or creates a new library if it do not exists.
    Returns the library id of the first library found with that name
    :param gi: Galaxy Instance
    :param library_name:
    :return: library Id
    '''
    libraries = library_list(gi)
    if library_name in [l['name'] for l in libraries]:
        return library_id(library_name)
    else:
        gi.libraries.create_library(library_name)
        return library_id(library_name)



def set_parameters(nstep, param_name, new_value):
    params = dict()
    params[nstep] = {param_name: new_value}
    return params
    #TODO delete


def set_params(json_wf, param_data):
    """
    Associate parameters to workflow steps via the step label. The result is a dictionary
    of parameters that can be passed to invoke_workflow.

    :param json_wf:
    :param param_data:
    :return:
    """
    params = {}
    for param_step_name in param_data:
        step_ids = (key for key, value in json_wf['steps'].items() if value['label'] == str(param_step_name))
        for step_id in step_ids:
            params.update({step_id: param_data[param_step_name]})
        for param_name in param_data[param_step_name]:
            if '|' in param_name:
                print("Workflow using Galaxy <repeat /> "
                                "param. type for {} / {}. "
                                "Make sure that workflow has as many entities of that repeat "
                                "as they are being set in the parameters file.".format(param_step_name, param_name))
                break
    return params


def load_input_files(gi, inputs, workflow, history):
    """
    Loads file in the inputs yaml to the Galaxy instance given. Returns
    datasets dictionary with names and histories. It associates existing datasets on Galaxy given by dataset_id
    to the input where they should be used.

    This setup currently doesn't support collections as inputs.

    Input yaml file should be formatted as:

    input_label_a:
      path: /path/to/file_a
      type:
    input_label_b:
      path: /path/to/file_b
      type:
    input_label_c:
      dataset_id:

    this makes it extensible to support
    :param gi: the galaxy instance (API object)
    :param inputs: dictionary of inputs as read from the inputs YAML file
    :param workflow: workflow object produced by gi.workflows.show_workflow
    :param history: the history object to where the files should be uploaded
    :return: inputs object for invoke_workflow
    """

    inputs_for_invoke = {}

    for step, step_data in workflow['inputs'].items():
        # upload file and record the identifier
        if step_data['label'] in inputs and 'path' in inputs[step_data['label']]:
            upload_res = gi.tools.upload_file(path=inputs[step_data['label']]['path'], history_id=history['id'],
                                              file_name=step_data['label'],
                                              file_type=inputs[step_data['label']]['type'])
            inputs_for_invoke[step] = {
                    'id': upload_res['outputs'][0]['id'],
                    'src': 'hda'
                }
        elif step_data['label'] in inputs and 'dataset_id' in inputs[step_data['label']]:
            inputs_for_invoke[step] = {
                'id': inputs[step_data['label']]['dataset_id'],
                'src': 'hda'
            }
        elif step_data['label'] in inputs and not isinstance(inputs[step_data['label']], Mapping):
            # We are in the presence of a simple parameter input
            inputs_for_invoke[step] = inputs[step_data['label']]
        else:
            raise ValueError("Label '{}' is not present in inputs yaml".format(step_data['label']))

    return inputs_for_invoke


# Nuevo
def set_params_json(json_wf, param_data):
    """
    Associate parameters to workflow steps via the step label. The result is a dictionary
    of parameters that can be passed to invoke_workflow.

    :param json_wf:
    :param param_data:
    :return:
    """
    params = {}
    for param_step_name in param_data:
        step_ids = (key for key, value in json_wf['steps'].items() if value['label'] == str(param_step_name))
        for step_id in step_ids:
            params.update({step_id: param_data[param_step_name]})
        for param_name in param_data[param_step_name]:
            if '|' in param_name:
                print("Workflow using Galaxy <repeat /> "
                                "param. type for {} / {}. "
                                "Make sure that workflow has as many entities of that repeat "
                                "as they are being set in the parameters file.".format(param_step_name, param_name))
                break
    return params


def validate_labels(wf_from_json, param_data, exit_on_error=True):
    """
    Checks that all workflow steps have labels (although if not the case, it will only
    warn that those steps won't be configurable for parameters) and that all step labels
    in the parameter files exist in the workflow file. If the second case is not true,
    the offending parameter is shown and the program exists.

    :param wf_from_json:
    :param param_data:
    :param exit_on_error:
    :return:
    """
    step_labels_wf = []
    for step_id, step_content in wf_from_json['steps'].items():
        if step_content['label'] is None:
            print("Step No {} in json workflow does not have a label, parameters are not mappable there.".format(step_id))
        step_labels_wf.append(step_content['label'])
    errors = 0
    for step_label_p, params in param_data.items():
        if step_label_p not in step_labels_wf:
            if exit_on_error:
                raise ValueError(
                    " '{}' parameter step label is not present in the workflow definition".format(step_label_p))
            errors += 1
    if errors == 0:
        print("Validation of labels: OK")


def validate_dataset_id_exists(gi, inputs):
    """
    Checks that dataset_id exists in the Galaxy instance when dataset_id are specified. Raises an error if the dataset
    id doesn't exists in the instance.

    :param gi:
    :param inputs:
    :return:
    """
    warned = False
    for input_key, input_content in inputs.items():
        if 'dataset_id' in input_content:
            ds_in_instance = gi.datasets.show_dataset(dataset_id=input_content['dataset_id'])
            if not warned:
                print("You are using direct dataset identifiers for inputs, "
                                "this execution is not portable accross instances.")
                warned = True
            if not isinstance(ds_in_instance, dict):
                raise ValueError("Input dataset_id {} does not exist in the Galaxy instance."
                                 .format(input_content['dataset_id']))


def validate_input_labels(wf_json, inputs):
    """
    Check that all input datasets in the workflow have labels, and that those labels are available in the inputs yaml.
    Raises an exception if those cases are not fulfilled.

    :param wf_json:
    :param inputs:
    :return: the number of input labels.
    """

    number_of_inputs = 0
    for step, step_content in wf_json['steps'].items():
        if step_content['type'] == 'data_input':
            number_of_inputs += 1
            if step_content['label'] is None:
                raise ValueError("Input step {} in workflow has no label set.".format(str(step)))

            if step_content['label'] not in inputs:
                raise ValueError("Input step {} label {} is not present in the inputs YAML file provided."
                                 .format(str(step), step_content['label']))
    return number_of_inputs


def validate_file_exists(inputs):
    """
    Checks that paths exists in the local file system.

    :param inputs: dictionary with inputs
    :return:
    """
    for input_key, input_content in inputs.items():
        if 'path' in input_content and not os.path.isfile(input_content['path']):
            raise ValueError("Input file {} does not exist for input label {}".format(input_content['path'], input_key))




def create_input(step: 'str', source:'str',id:'str'):
    '''
    Create a workflow input
    :param step:  step input
    :param source: 'hdda' (history); 'lbda' (library)
    :param id: dataset id
    :return: Input dictionary
    '''
    datamap = dict()
    datamap[step] = {'id': id, 'src': source}
    return datamap



def run_workflow(gi , name: 'str', input_path: 'str',input_name, *, step_index: 'str' ='0', history_name: 'str' ='Test_History', params = None)->dict:
    '''
    Directly run a workflow from a local file
    #1 Create New History
    #2 Upload Dataset
    #3 create the new input using upload_file tool in the History just created
    #4 Invoke the Workflow in the new History -> this Generates: id (as WorkflowInvocation), workflow_id, history_id,
                                                and a Job id for every step executed.

    :param gi: galaxy Instance
    :param name: Workflow Name
    :param input_path: Input path
    :param input_name: Name
    :param step_index: Input step index
    :param history_name: Name of the history where the workflow will be invoked. If this history does not exist
    it will be created.
    :return:an invocation dictionary
    '''
    h_id = gi.histories.create_history(name=history_name)['id']
    d_id = gi.tools.upload_file(input_path, h_id, filename= input_name)['outputs'][0]['id']
    w_id = workflow_id(gi, name)
    dataset = {'src': 'hda', 'id': d_id}
    invocation = gi.workflows.invoke_workflow(w_id,
                                              inputs ={step_index: dataset},
                                              history_id= h_id,
                                              inputs_by ='step_index',
                                              params = params
                                              )
    return invocation

def run_workflow_files(gi,wf_name,input_file: 'str',param_file: 'str',history_name):
    param_data = read_yaml_file(param_file)
    inputs_data = read_yaml_file(input_file)
    w_id = workflow_id(gi, wf_name)
    wf_dict = gi.workflows.export_workflow_dict(w_id)
    show_wf = gi.workflows.show_workflow(w_id)
    allowed_error_states = {'tools': {}, 'datasets': set()}
    # Move any simple parameters from parameters to inputs
    params_to_move = []
    for pk, pv in param_data.items():
        if not isinstance(pv, dict):
            #TODO params could be nested dictionary?
            params_to_move.append(pk)

    for pk in params_to_move:
        inputs_data[pk] = param_data[pk]

    validate_labels(wf_dict, param_data)
    num_inputs = validate_input_labels(wf_json=wf_dict, inputs=inputs_data)
    if num_inputs > 0:
        validate_file_exists(inputs_data)

    validate_dataset_id_exists(gi, inputs_data)

    print('Create new history to run workflow ...')
    if num_inputs > 0:
        history = gi.histories.create_history(name=history_name)
        datamap = load_input_files(gi, inputs=inputs_data,
                                   workflow=show_wf, history=history)
        print('Set parameters ...')
        params = set_params(wf_dict, param_data)

    print('Running workflow {}...'.format(show_wf['name']))
    results = gi.workflows.invoke_workflow(workflow_id=w_id,
                                           inputs=datamap,
                                           params=params,
                                           history_name=history_name)

    # time.sleep(100)
    return results #invocation objectt

def get_job(gi,invocation,step):
    '''
    Job information from an invocation given the step of interest. A job is the execution of a step (in a workflow)
    :param gi:
    :param invocation:
    :return:
    '''
    step = gi.invocations.show_invocation(invocation['id'])['steps'][int(step)]
    state = step['state']
    job_id = step['job_id']
    job =  gi.jobs.show_job(job_id)
    return job


def invocation_percent_complete(gi,invocation)->'int':
    status = gi.histories.get_status(invocation['history_id'])
    return status['percent_complete']
# TODO WRITE A BETER WAY TO CHECK PROCESSES

def invocation_errors(gi,invocation)->'int':
    status = gi.histories.get_status(invocation['history_id'])
    return status['state_details']['error']



def list_invocation_results(gi,invocation_id: 'str'):
    '''
    Generates a list of results from an invocation
    :param gi: Galaxy Instance
    :param invocation_id: Workflow Invocation Id
    :return: a list of dictionary of datasets or a error
    '''
    state = gi.invocations.show_invocation(invocation_id)['state']
    if state == 'scheduled':
        h_id = gi.invocations.show_invocation(invocation_id)['history_id']
        if gi.histories.get_status(h_id)['state'] == 'ok':
            results = gi.histories.show_matching_datasets(h_id)
            return results
        else:
            return gi.histories.get_status(h_id)




def download_result(gi, results:'list', path:'str'):
    '''
    Download invocation result
    :param gi:
    :param results:
    :param path:
    :return:
    '''
    if isinstance(results,list):
        for r in results:
            gi.datasets.download_dataset(r['id'], file_path=path)
    else:
        print(results)
    # TODO delete history after successfully download


def export_workflow(instance, wf_name):
    gi = instance.connect()
    wf_id = get_workflow_from_name(gi,wf_name)
    wf_dic = gi.workflows.import_workflow_dict(wf_id)
    return wf_dic

def import_workflow(instance, wf_Id):
    gi = instance.connect()
    wf = gi.workflows.import_workflow_dict(wf_Id)
    return wf


def check_tools(wf1_dic,wf2_dic):
    steps1 = wf1_dic['step']
    steps2 = wf2_dic['step']
    tool_list = list()
    for step, content in steps1.items():
        if step['errors'] == "Tool is not installed":
            tool_list.append(steps2[step]['tool_shed_repository'])
    if len(tool_list) == 0:
        return 'all tools are instaled'
    else:
        return tool_list


def install_tools(instance,tools):

    '''
    tool_shed_url, name, owner,
    changeset_revision,
    install_tool_dependencies = False,
    install_repository_dependencies = False,
    install_resolver_dependencies = False,
    tool_panel_section_id = None,
    new_tool_panel_section_label = Non
    '''
    if isinstance(tools, list):
        gi = instance.connect()
        for tool in tools:
            tool_shed_url = 'https://' + tool['tool_shed']
            gi.tooshed.install_repository_revision(tool_shed_url=tool_shed_url,
                                                   name=tool['name'],
                                                   owner=tool['owner'],
                                                   changeset_revision= tool['changeset_revision'],
                                                   install_tool_dependencies=True,
                                                   install_repository_dependencies=True,
                                                   install_resolver_dependencies=True,
                                                   new_tool_panel_section_label='New'
                                                   )

class instance():
    def __init__(self,name,file):
        self.name = name
        self.file = file


    def connect(self):
        gi = galaxy_instance(self.file, name=self.name)
        return gi














