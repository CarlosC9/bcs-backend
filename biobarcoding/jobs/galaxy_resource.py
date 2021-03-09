from typing import List, Any

from bioblend import galaxy
import json
import yaml
import os
import re

from biobarcoding.common import ROOT
from biobarcoding.db_models import DBSession
from biobarcoding.db_models.jobs import ComputeResource
from biobarcoding.jobs import JobExecutorAtResource
from pathlib import Path
from multiprocessing import Process, Queue



# no tiene sentido que haga esto xq galaxy ya es este objeto
class instance():
    def __init__(self, name, file):
        self.name = name
        self.file = file

    def connect(self):  # esta función es la que debería crearme la instance
        gi = galaxy_instance(self.file, name=self.name)
        return gi


def galaxy_instance(path, name='__default'):
    data = read_yaml_file(path)
    assert name in data, 'unknown instance'
    gal = data[name]
    if isinstance(gal, dict):
        return gal
    else:
        return data[gal]


def login(api_key: 'str', url: 'str') -> object:
    user_key = api_key
    gi = galaxy.GalaxyInstance(url=url, key=user_key)
    return gi


def library_list(gi):
    libraries = gi.libraries.get_libraries()
    return libraries


def library_id(gi, libname: 'str'):
    libraries = gi.libraries.get_libraries()
    try:
        lib = next(item for item in libraries if item["name"] == libname)
    except StopIteration:
        return 'there is no library named{}'.format(libname)
    else:
        return lib['id']


def get_history_id(gi, history_name: 'str'):
    histories = gi.histories.get_histories()
    try:
        hist = next(item for item in histories if item["name"] == history_name)
    except StopIteration:
        return 'there is no library named{}'.format(history_name)
    else:
        return hist['id']


def workflow_list(gi):
    workflows = gi.workflows.get_workflows()
    return workflows


def workflow_id(gi, name: 'str'):
    workflows = gi.workflows.get_workflows()
    try:
        work = next(item for item in workflows if item["name"] == name)
    except StopIteration:
        return 'there is no workflow named{}'.format(name)
    else:

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


def parse_input(input: 'str', ext: 'str'):
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

def delete_history(gi, history):
    gi.histories.delete_histoy(history['id'])


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
        for input in inputs:
            if step_data['label'] in input.values() and 'path' in input.keys():
                d_id = gi.histories.show_matching_datasets(history_id=history['id'],name_filter=step_data['label'])[0].get('id')
                inputs_for_invoke[step] = {
                    'id': d_id,
                    'src': 'hda'
                }
    return inputs_for_invoke


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
            print(
                "Step No {} in json workflow does not have a label, parameters are not mappable there.".format(step_id))
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


def create_input(step: 'str', source: 'str', id: 'str'):
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


def params_input_creation(gi, workflow_name, inputs_data, param_data, history_id=None, history_name=None):
    '''
        set and chel params dictionaty and load data and set datamap

    '''
    w_id = workflow_id(gi, workflow_name)
    wf_dict = gi.workflows.export_workflow_dict(w_id)
    show_wf = gi.workflows.show_workflow(w_id)
    params_to_move = []
    for pk, pv in param_data.items():
        if not isinstance(pv, dict):
            # TODO params could be nested dictionary?
            params_to_move.append(pk)

    for pk in params_to_move:
        inputs_data[pk] = param_data[pk]
    history = gi.histories.get_histories(history_id)[0]
    datamap = load_input_files(gi, inputs=inputs_data,
                               workflow=show_wf, history=history)
    # TODO check that input parameters are correct by form
    print('Set parameters ...')
    params = set_params(wf_dict, param_data)
    return datamap, params


def params_creation(gi, workflow_name, param_data):  # TODO test
    w_id = workflow_id(gi, workflow_name)
    wf_dict = gi.workflows.export_workflow_dict(w_id)
    params_to_move = []
    for pk, pv in param_data.items():
        if not isinstance(pv, dict):
            params_to_move.append(pk)

    validate_labels(wf_dict, param_data)
    params = set_params(wf_dict, param_data)
    return params


def input_creation(gi, workflow_name, history_name, inputs_data):  # TODO test
    w_id = workflow_id(gi, workflow_name)
    wf_dict = gi.workflows.export_workflow_dict(w_id)
    show_wf = gi.workflows.show_workflow(w_id)
    num_inputs = validate_input_labels(wf_json=wf_dict, inputs=inputs_data)
    if num_inputs > 0:
        validate_file_exists(inputs_data)

    validate_dataset_id_exists(gi, inputs_data)

    print('Create new history to run workflow ...')
    history = gi.histories.create_history(name=history_name)
    #dict step : data id
    datamap = load_input_files(gi, inputs=inputs_data,
                               workflow=show_wf, history=history)
    # TODO check that input parameters are correct
    return datamap


def run_workflow_files(gi, wf_name, input_file: 'str', param_file: 'str', history_name):
    param_data = read_yaml_file(param_file)
    inputs_data = read_yaml_file(input_file)
    w_id = workflow_id(gi, wf_name)
    show_wf = gi.workflows.show_workflow(w_id)
    datamap, parameters = params_creation(gi, wf_name, inputs_data, param_data)
    print('Running workflow {}...'.format(show_wf['name']))
    invocation = gi.workflows.invoke_workflow(workflow_id=w_id,
                                              inputs=datamap,
                                              params=parameters,
                                              history_name=history_name)
    return invocation


def get_historyID_by_invocation(gi, h_id):
    invocations = gi.invocations.get_invocations()
    try:
        inv = next(item for item in invocations if item["history_id"] == h_id)
    except StopIteration:
        return 'there is no invocation named{}'.format(h_id)
    else:
        return inv['id']


def invocation_errors(gi, invocation_id) -> 'dict':
    invocation = gi.invocations.show_invocation(invocation_id)
    status = gi.histories.get_status(invocation['history_id'])
    state = status['state']
    if state == 'error':
        error_detail = status['state_details']
        return error_detail
    else:
        return state


def get_job(gi, invocation, step):
    '''
    Job information from an invocation given the step of interest. A job is the execution of a step (in a workflow)
    :param gi:
    :param invocation:
    :return:
    '''
    step = gi.invocations.show_invocation(invocation['id'])['steps'][int(step)]
    state = step['state']
    job_id = step['job_id']
    job = gi.jobs.show_job(job_id)
    return job


def list_invocation_results(gi, invocation_id: 'str'):
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


def download_result(gi, results: 'list', path: 'str'):
    '''
    Download invocation result
    :param gi:
    :param results:
    :param path:
    :return:
    '''
    if isinstance(results, list):
        for r in results:
            gi.datasets.download_dataset(r['id'], file_path=path)
    else:
        print(results)
    # TODO delete history after successfully download


def export_workflow(instance, wf_name):
    gi = instance.connect()
    wf_id = get_workflow_from_name(gi, wf_name)
    wf_dic = gi.workflows.import_workflow_dict(wf_id)
    return wf_dic


def import_workflow(instance, wf_Id):
    gi = instance.connect()
    wf = gi.workflows.import_workflow_dict(wf_Id)
    return wf


def check_tools(wf1_dic, wf2_dic):
    steps1 = wf1_dic['steps']
    steps2 = wf2_dic['steps']
    tool_list = list()
    for step, content in steps1.items():
        if 'errors' in content:
            if content['errors'] == "Tool is not installed":
                # TODO depende de la versión de galaxi esto lleva un punto al final o no xq lo que hay que buscar
                #  otra cosa
                tool_list.append(steps2[step]['tool_shed_repository'])
    if len(tool_list) == 0:
        return 'all tools are installed'
    else:
        return tool_list


def install_tools(gi, tools):
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
        for tool in tools:
            tool_shed_url = 'https://' + tool['tool_shed']
            gi.toolshed.install_repository_revision(tool_shed_url=tool_shed_url,
                                                    name=tool['name'],
                                                    owner=tool['owner'],
                                                    changeset_revision=tool['changeset_revision'],
                                                    install_tool_dependencies=True,
                                                    install_repository_dependencies=True,
                                                    install_resolver_dependencies=True,
                                                    new_tool_panel_section_label='New'
                                                    )


class JobExecutorAtGalaxy(JobExecutorAtResource):
    def __init__(self):
        self.api_key = None
        self.url = None
        self.galaxy_instance = None

    def set_resource(self, params):
        self.api_key = params['jm_credentials']['api_key']
        self.url = params['jm_location']['url']

    def connect(self):
        self.galaxy_instance = login(self.api_key, self.url)

    def disconnect(self):
        # todo
        pass

    def create_job_workspace(self, name):
        self.connect()
        gi = self.galaxy_instance
        history = gi.histories.create_history(name=str(name))
        # tengo que retornar algo diferente si no se puede conectar a galaxy
        self.disconnect()
        return history['id']

    def remove_job_workspace(self, workspace):
        self.connect()
        gi = self.galaxy_instance
        # TODO hacer un purge y tmb haer un purge del dataset..... p quizás poner como tarea de mantenimiento del celery??
        gi.histories.delete_history(get_history_id(gi, workspace))

    def upload_file(self,local_path, **kwards):
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

        self.connect()
        gi = self.galaxy_instance
        inputs = 'files list'
        workflow = gi.workflows.show_workflow
        inputs_for_invoke = {}
        workflow_id = kwards.get('workflow')
        history = kwards.get('workspace')
        h_id = get_history_id(gi,history)
        label = kwards.get('step')
        try:
            upload_info = gi.tools.upload_file(path=local_path,
                                               history_id=h_id,
                                               file_name=label)
            pid = upload_info['jobs'][0]['id']
        except:
            pid = "error"
        return pid


    def exists(self, **kwargs):
        self.connect()
        gi = self.galaxy_instance
        label = kwargs.get('step')
        history = kwargs.get('workspace')
        h_id = get_history_id(gi,history_name= history)
        dataset_info = gi.histories.show_matching_datasets(history_id = h_id , name_filter = label )
        if len(dataset_info)>0:
            return True
        else:
            return False



    def submit(self, workspace, params):
        # "process": {"name": "MSA ClustalW",
        #             "inputs":
        #                 {"parameters":
        #                      {"MSA ClustalW": {"darna": "PROTEIN"}
        #                       },
        #                  "data": [{"step":"Input dataset" ,
        #                           "path": "/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/matK_25taxones_Netgendem_SINalinear.fasta",
        #                           "type": "fasta"},....]
        #                  }
        #             }

        self.connect()
        gi = self.galaxy_instance
        input_params = params['inputs']['parameters']
        inputs = params['inputs']['data']
        workflow = params['name']
        w_id = workflow_id(gi, workflow)
        datamap, parameters = params_input_creation(gi, workflow, inputs, input_params,
                                                    history_name=workspace)  # catch error
        # TODO revisar la fomra del data map, quizás me pueda ahorrar toda esa función xq esto ya está comprobadoo desde el gui
        history_id = gi.histories.get_histories(name=workspace)[0]['id']
        invocation = gi.workflows.invoke_workflow(workflow_id=w_id,
                                                  inputs=datamap,
                                                  params=parameters,
                                                  history_id=history_id)
        return invocation['id']

    def job_status(self, pid):
        """
        :return: state of the given job among the following values: `new`,
          `queued`, `running`, `waiting`, `ok`. If the state cannot be
          retrieved, an empty string is returned.
        """
        self.connect()
        gi = self.galaxy_instance
        if pid:
            state = gi.jobs.get_state(pid)
            # TODO al pasarle el id del invocation estoy viendo el estado de todo el workflow? una invocación de
            #  workflow debe ser un job a su vez.....
            print(f"{state} job in job_status function")
            return state
        else:
            return None

    def cancel_job(self, native_id):
        self.connect()
        gi = self.galaxy_instance
        gi.invocations.cancel_invocation(native_id)
        # job here refers to invocation

    def download_file(self, native_id,local_path = None):
        self.connect()
        gi = self.galaxy_instance
        # todo bajar documento uno a uno
        # 1. cargar lista
        # 2. comprobar que el documento ya esté
        # 3. si no está descargarlo al path
        r = list_invocation_results(gi, native_id)
        download_result(gi, r, '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/download')
        return r

def convert_workflows_to_formly():
    wfdict1 = {'clustalw': ROOT + '/biobarcoding/inputs_schema/clustalw_galaxy.json',
               'phyml': ROOT + '/biobarcoding/inputs_schema/phyml_galaxy.json',
               'fname' : 'clustalw_phyml_formly.json'
               }
    wfdict2 = {'clustalw': ROOT + '/biobarcoding/inputs_schema/clustalw_galaxy.json',
               'fname' : 'clustalw_formly.json'}
    path = ROOT + '/biobarcoding/inputs_schema/'
    lwdict = [wfdict1, wfdict2]
    for wfdict in lwdict:
        if wfdict['fname'] not in os.listdir(path):
            newpath = path + wfdict['fname']
            del wfdict['fname']
            convertToFormly(wfdict,newpath)

def initialize_galaxy(flask_app):
    if {'GALAXY_API_KEY', 'GALAXY_LOCATION'} <= flask_app.config.keys():
        api_key = flask_app.config['GALAXY_API_KEY']
        url = flask_app.config['GALAXY_LOCATION']

        # Update resource location
        session = DBSession()
        local_uuid = "8fac3ce8-8796-445f-ac27-4baedadeff3b"
        r = session.query(ComputeResource).filter(ComputeResource.uuid == local_uuid).first()
        if r:
            r.jm_location = {"url": url}
            r.jm_credentials = {"api_key": api_key}
            session.commit()
        DBSession.remove()

        # Install basic workflow if it is not installed
        gi = login(api_key, url)
        path = ROOT + '/biobarcoding/workflows/'
        for workflow in os.listdir(path):
            workflow_path = path + workflow
            with open(workflow_path, 'r') as f:
                wf_dict_in = json.load(f)
            name = wf_dict_in['name']
            if workflow_id(gi, name) == 'there is no workflow named{}'.format(name):
                wf = gi.workflows.import_workflow_from_local_path(workflow_path)
                wf_dict_out = gi.workflows.export_workflow_dict(wf['id'])
                list_of_tools = check_tools(wf_dict_out, wf_dict_in)
                install_tools(gi, list_of_tools)
        # conversion of workflows galaxy into workflows formly
        convert_workflows_to_formly()
    else:
        return 'No Galaxy test credentials in config file'


class ToFormlyConverter:
    field = {
        # 'original' : 'formly'
        'name': 'key',
        'type': 'type'
    }
    templateoptionsfields = {
        # 'original' : 'formly'
        'label': 'label',
        'optional': 'required'
    }

    def __init__(self, step_label):
        self.step_label = step_label

    @staticmethod
    def rename_keys(d, keys):
        return dict([(keys.get(k), v) for k, v in d.items() if k in keys.keys()])

    @staticmethod
    def options(g_input):
        return [{'value': o[1], 'label': o[0]} for o in g_input['options']]

    def choose_converter(self, g_input):
        input_type = g_input['model_class']
        if input_type == 'SelectToolParameter':
            converter = convertSelectToolParameter(self.step_label)
        elif input_type == 'BooleanToolParameter':
            converter = convertBooleanToolParameter(self.step_label)
        elif input_type == 'Conditional':
            converter = converterConditional(self.step_label)
        elif input_type == 'IntegerToolParameter':
            converter = converterIntegerToolParameter(self.step_label)
        elif input_type == 'FloatToolParameter':
            converter = converterIntegerToolParameter(self.step_label)
        elif input_type == 'TextToolParameter':
            converter = converterTextToolParameter(self.step_label)
        else:
            return 'no converter for model class {}'.format(g_input['model_class'])
        return converter.convert(g_input)

    def conversion(self, g_input):
        form = self.rename_keys(g_input, self.field)
        form['key'] = '.'.join([self.step_label, form['key']])
        form['templateOptions'] = self.rename_keys(g_input, self.templateoptionsfields)
        if 'value' in g_input:
            form['defaultValue'] = g_input['value']
        form['templateOptions']['required'] = not form['templateOptions']['required']
        form['templateOptions']['description'] = g_input['help']
        return form

    def get_formly_dict(self, g_input):
        l = []
        for i in g_input:
            forms = self.choose_converter(i)
            if isinstance(forms, list) and len(forms) > 1:
                for f in forms:
                    if isinstance(f, dict):
                        l.append(f)
            else:
                if isinstance(forms, dict):
                    l.append(forms)
        return l


class convertBooleanToolParameter(ToFormlyConverter): #TODO hay algún fallo al poner el valor por defecto!
    def __init__(self,step_label):
        super(convertBooleanToolParameter, self).__init__(step_label)

    def convert(self, g_input):
        form = self.conversion(g_input)
        form['type'] = 'radio'
        form['templateOptions']['options'] = [
            {'value': g_input['truevalue'], 'label': 'Yes'},
            {'value': g_input['falsevalue'], 'label': 'No'}  # check is needed
        ]
        if form['defaultValue']:
            if form['defaultValue'] == 'false':
                form['defaultValue'] = 'OFF'
            elif form['defaultValue'] == 'true':
                form['defaultValue'] = 'YES'
            else:
                form['defaultValue'] = form['defaultValue']

        return form


class convertSelectToolParameter(ToFormlyConverter):
    def __init__(self, step_label):
        super(convertSelectToolParameter, self).__init__(step_label)

    def convert(self, g_input):
        form = self.conversion(g_input)
        form['templateOptions']['options'] = self.options(g_input)
        return form

class converterIntegerToolParameter(ToFormlyConverter):
    def __init__(self,step_label):
        super(converterIntegerToolParameter, self).__init__(step_label)

    def convert(self, g_input):
        form = self.conversion(g_input)
        form['type'] = 'input'
        form['templateOptions']['type'] = 'number'
        if g_input['min'] != None:
            form['templateOptions']['min'] = int(g_input['min'])
        if g_input['max'] != None:
            form['templateOptions']['max'] = int(g_input['max'])
        return form


class converterTextToolParameter(ToFormlyConverter):
    def __init__(self,step_label):
        super(converterTextToolParameter, self).__init__(step_label)
    def convert(self, input):
        form = self.conversion(input)
        form['type'] = 'textarea'
        return form


class converterConditional(ToFormlyConverter):
    def __init__(self,step_label):
        super(converterConditional, self).__init__(step_label)
        self.selector = None
        self.cases = None

    def convert(self, g_inputs):
        self.selector = g_inputs['test_param']
        self.cases = g_inputs['cases']
        form = list()
        selector_form = self.choose_converter(self.selector)
        form.append(selector_form)
        for i in self.cases:
            if len(i['inputs']) > 0:
                for j in i['inputs']:
                    case_form = self.choose_converter(j)
                    if isinstance(case_form, dict):
                        case_form['hideExpression'] = 'model.' + selector_form['key'] + '!=\'' + i['value'] + '\''
                    else:
                        print('no converter for ', j['model_class'])
                    form.append(case_form)
        return form


def convertToFormly(wf_steps, newpath):
    '''
    dicctionary step_label: form_path
    '''
    fieldGroup = list()
    formly = dict()
    formly['type'] = 'stepper'
    for k, v in wf_steps.items():
        input_path = v
        with open(input_path, 'r') as f:
            galaxy_dict_in = json.load(f)
        inputs = galaxy_dict_in['inputs']
        convert = ToFormlyConverter(k)
        form = dict()
        form['templateOptions'] = {'label': k}
        form['fieldGroup'] = convert.get_formly_dict(inputs)
        fieldGroup.append(form)
    formly['fieldGroup'] = fieldGroup
    formly_json = json.dumps([formly], indent=3)
    with open(newpath, 'w') as file:
        file.write(formly_json)
