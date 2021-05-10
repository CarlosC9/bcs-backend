import time
from urllib.parse import urljoin

from bioblend import galaxy
import json

from biobarcoding.common import ROOT
from biobarcoding.db_models import DBSession
from biobarcoding.db_models.jobs import ComputeResource
from biobarcoding.jobs import JobExecutorAtResource

import os

# GALAXY HELPERS
from biobarcoding.tasks.definitions import write_to_file


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
        return None
    else:
        return hist['id']


def workflow_list(gi):
    workflows = gi.workflows.get_workflows()
    return workflows


def workflow_id(gi, name: 'str'):
    workflows = gi.workflows.get_workflows()
    try:
        workflow = next(item for item in workflows if item["name"] == name)
    except StopIteration:
        return None
    else:

        return workflow['id']


def import_workflow_from_file(gi, workflow_file):
    imported_workflow = [gi.workflows.import_workflow_from_local_path(file_local_path=workflow_file)]
    return imported_workflow


def get_workflow_from_name(gi, workflow_name):
    workflow = gi.workflows.get_workflows(name=workflow_name)
    return workflow


def get_workflow_id(wf):
    for wf_dic in wf:
        wf_id = wf_dic['id']
    return wf_id


def dataset_list(gi):
    datasets = gi.datasets.get_datasets()
    return datasets


def dataset_id(gi, name: 'str'):
    datasets = dataset_list(gi)
    dataset = next(item for item in datasets if item["name"] == name)
    return dataset['id']


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
        for step_id in list(step_ids):
            params.update({step_id: param_data[param_step_name]})
        for param_name in param_data[param_step_name]:
            if '|' in param_name:
                print("Workflow using Galaxy <repeat /> "
                      "param. type for {} / {}. "
                      "Make sure that workflow has as many entities of that repeat "
                      "as they are being set in the parameters file.".format(param_step_name, param_name))
                break
    return params


def get_datamap(gi, inputs, workflow, history):
    '''
    Associate input labels with dataset id that already exists in a history
    @param gi: galaxy instance
    @param inputs: list of inputs labels
    @param workflow: workflow dictionary
    @param history: history dictionary
    @return:  {  id: galaxy dataset id,
                'src:'hda } * That means that the dataset is stored in a history.
    '''

    inputs_for_invoke = {}

    for step, step_data in workflow['inputs'].items():
        for input in inputs:
            if step_data['label'] in input.values() and 'remote_name' in input.keys():
                d_id = gi.histories.show_matching_datasets(history_id=history['id'],name_filter=step_data['label'])[0].get('id')
                inputs_for_invoke[step] = {
                    'id': d_id,
                    'src': 'hda'
                }
    return inputs_for_invoke


def params_input_creation(gi, workflow_name, inputs_data, param_data, history_id=None, history_name=None):
    # TODO SEPARAR LA FUNCIÓN DE GET_DATAMAP Y PASAR AL ADAPTOR
    '''
        set and check params dictionary and load data and set datamap

    '''
    print(f"inputs_data: {inputs_data}")
    print(f"Param data: {param_data}")
    w_id = workflow_id(gi, workflow_name)
    print(f"Workflow id: {w_id}")
    wf_dict = gi.workflows.export_workflow_dict(w_id)
    show_wf = gi.workflows.show_workflow(w_id)
    print(f"Show Workflow: {show_wf}")
    params_to_move = []
    for pk, pv in param_data.items():
        if not isinstance(pv, dict):
            # TODO params could be nested dictionary?
            params_to_move.append(pk)

    for pk in params_to_move:
        inputs_data[pk] = param_data[pk]
    history = gi.histories.get_histories(history_id)[0]
    datamap = get_datamap(gi, inputs=inputs_data, workflow=show_wf, history=history)
    print('Set parameters ...')
    params = set_params(wf_dict, param_data)
    return datamap, params

def get_history_id_by_invocation(gi, h_id):
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


def get_job_from_invocation(gi, invocation, step):
    '''
    Job information from an invocation given the step of interest. A job is the execution of a step (in a workflow)
    An invocation is the result of calling a pre-defined workflow. So an invocation is composed of several jobs.
    :param : galaxy instance
    :param invocation: invocation dictionary.
    :return: A job dictionary
    '''
    step = gi.invocations.show_invocation(invocation['id'])['steps'][int(step)]
    #TODO repasar este metodo
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



def get_stdout_stderr(gi,result,history_name):
    remote_name = result['remote_name']
    dataset_id = \
    gi.histories.show_matching_datasets(history_id=get_history_id(gi, history_name), name_filter=remote_name)[
        0]['id']
    provenance = gi.histories.show_dataset_provenance(history_id=get_history_id(gi, history_name),
                                                      dataset_id=dataset_id)
    std = dict(stderr = provenance['stderr'],stdout = provenance['stdout'])
    return std



class JobExecutorAtGalaxy(JobExecutorAtResource):
    def __init__(self, job_id):
        super().__init__(job_id)
        self.api_key = None
        self.url = None
        self.galaxy_instance = None
        self.workspace = job_id

    def set_resource(self, params):
        self.api_key = params['jm_credentials']['api_key']
        self.url = params['jm_location']['url']


    def connect(self):
        self.galaxy_instance = login(self.api_key, self.url)

    def check(self):
        self.connect()
        gi = self.galaxy_instance
        try:
            gi.config.get_config()
        except: # ConnectionError dosnt work
            return False
        return True

    def disconnect(self):
        # necessary?
        pass

    def get_upload_files_list(self, job_context):
        return job_context["process"]["inputs"]["data"]

    def get_download_files_list(self, job_context):
        return job_context["results"]

    def create_job_workspace(self, name):
        self.connect()
        gi = self.galaxy_instance
        history = gi.histories.create_history(name=str(name))
        # tengo que retornar algo diferente si no se puede conectar a galaxy
        self.disconnect()
        return history['id']

    def job_workspace_exists(self, job_id):
        self.connect()
        gi = self.galaxy_instance
        history_id = get_history_id(gi,str(job_id))
        if history_id:
            return True
        else:
            return False

    def remove_job_workspace(self, job_context):
        self.connect()
        gi = self.galaxy_instance
        history_name = self.workspace
        gi.histories.delete_history(get_history_id(gi, history_name),purge = True)

    def upload_file(self,job_context):
        self.connect()
        gi = self.galaxy_instance
        i = job_context["state_dict"]["idx"]
        label = job_context["process"]["inputs"]["data"][i]["remote_name"]
        local_path = job_context["process"]["inputs"]["data"][i]["file"]
        history = self.workspace
        upload_info = gi.tools.upload_file(path=local_path,
                                           history_id=get_history_id(gi,history),
                                           file_name=label)
        pid = upload_info['jobs'][0]['id']
        return pid

    def exists(self,job_context):
        self.connect()
        gi = self.galaxy_instance
        if not job_context.get('state_dict'):
            return False
        i = job_context["state_dict"]["idx"]
        if job_context["state_dict"]["state"] == "upload":
            # check if exists locally
            file_local_path = os.path.join(self.LOCAL_WORKSPACE, self.workspace,
                                           self.get_upload_files_list(job_context)[i]["file"])
            try:
                local_size = os.path.getsize(file_local_path)
            except FileNotFoundError:
                print(f"File {file_local_path} not found in your local system")
                return None
            # check if exists remotely
            label = self.get_upload_files_list(job_context)[i]["remote_name"]
            h_id = get_history_id(gi,history_name= self.workspace)
            dataset_info = gi.histories.show_matching_datasets(history_id = h_id , name_filter = label )
            # check that are the same
            if len(dataset_info) > 0:
                self.__write_logs(job_context)
                remote_size = dataset_info[0]['file_size']
                # TODO not implemented: there is a little difference between sizes
                # if local_size == remote_size:
                return True
        if job_context["state_dict"]["state"] == "download":
            file_local_path = os.path.join(self.LOCAL_WORKSPACE,self.workspace,self.get_download_files_list(job_context)[i]["file"])
            if os.path.exists(os.path.join(file_local_path)):
                return True
            else:
                print(f"File {file_local_path} not found in your local system")
                return None



    def submit(self, job_context):
        self.connect()
        gi = self.galaxy_instance
        input_params = job_context['inputs']['parameters']
        inputs = job_context['inputs']['data']
        workflow = job_context['workflow_name']
        w_id = workflow_id(gi, workflow)
        workspace = self.workspace
        datamap, parameters = params_input_creation(gi, workflow, inputs, input_params,
                                                    history_name=workspace)
        history_id = gi.histories.get_histories(name=workspace)[0]['id']
        invocation = gi.workflows.invoke_workflow(workflow_id=w_id,
                                                  inputs=datamap,
                                                  params=parameters,
                                                  history_id=history_id)
        print(f"Invocation: {invocation}")
        return invocation['id']

    def __invocation_status(self, job_context):
        """
        Returns the state of this history

        :type history_id: str
        :param history_id: Encoded history ID

        :rtype: dict
        :return: A dict documenting the current state of the history. Has the following keys:
            'state' = This is the current state of the history, such as ok, error, new etc.
            'state_details' = Contains individual statistics for various dataset states.
            'percent_complete' = The overall number of datasets processed to completion.
        """
        gi = self.galaxy_instance
        status = gi.histories.get_status(get_history_id(gi,str(self.workspace)))
        write_to_file(self.log_filenames_dict['submit_stdout'],json.dumps(status['state_details']))
        write_to_file(self.log_filenames_dict['submit_stdout'], 'percent complete: ' + json.dumps(status['percent_complete']))
        "the first call to state cab be ok just because is the state of previous job. "
        print(f"{status} job in job_status function")
        files_list = job_context['process']['inputs']['data']
        if status['state'] == 'ok' and status['state_details']['ok'] == len(files_list):
            # galaxy is still displaying status for the upload process
            status['state'] = "queued"
            print(f"changing PROCESS status from ok to queued this should happend only a few times")
        elif status['state'] == 'error':
            invocation = self.galaxy_instance.invocations.show_invocation(job_context["pid"])
            for step in invocation["steps"]:
                if step["job_id"] is not None:
                    job_context["pid"] = step["job_id"]
                    self.__write_logs(job_context)
            return status['state_details']
        elif status['state'] == 'ok':
            invocation = self.galaxy_instance.invocations.show_invocation(job_context["pid"])
            for step in invocation["steps"]:
                if step["job_id"] is not None:
                    job_context["pid"] = step["job_id"]
                    self.__write_logs(job_context)
        return status['state']

    def __upload_status(self,job_context):
        gi = self.galaxy_instance
        i = job_context["state_dict"]["idx"]
        if i == 0:
            # todavía no he empezado
            return None
        status = gi.histories.get_status(get_history_id(gi, job_context['pid']))
        if status['state'] == 'ok' and status['state_details']['ok'] < len(i):
            # "in this case we can say that galaxy does not know about de nwe job"
            status['state'] = "queued"
            print(f"changing {i} UPLOAD status from ok to queued this should happend only a few times")

        if status['state'] == 'error':
            print('writing log....')
            self.__write_logs(job_context)  # me da que puede ser lo mismo
            return status['state_details']
        else:
            print('writing log....')
            self.__write_logs(job_context)
            return status['state']


    def step_status(self, job_context):
        """
           :return: state of the given job among the following values: `new`,
             `queued`, `running`, `waiting`, `ok`. If the state cannot be
             retrieved, an empty string is returned.
           """
        self.connect()
        if job_context.get("state_dict"):
            if job_context["state_dict"].get("state") == "upload":
                status = self.__upload_status(job_context) # i need
                return status # if error return dict
            if job_context["state_dict"].get("state") =="download":
                pid = job_context["pid"]
                return self.local_job_status(self.workspace,pid)
            if job_context["state_dict"].get("state") == "submit":
                status = self.__invocation_status(job_context)
                return status
        else:
            return None


    def cancel_job(self, native_id):
        self.connect()
        gi = self.galaxy_instance
        gi.invocations.cancel_invocation(native_id)
        # job here refers to invocation


    def __write_logs(self,job_context):
        '''
        write stdout and stderr  in <state>_staout and stderr files:
        when state = upload stdout will be always empty
        when state = submit stdout is retrieved several times to check whether there is a latency to retriece std
        from galaxy.
        output dictionary when asking for upload tool job
        'outputs': {'
            {'output0':
                {'id': 'f2db41e1fa331b3e',
               'src': 'hda',
               'uuid': '789cf547-a1a6-489d-8cb5-ab159ba9c384'}}

       Output dictionary when job id refers to an invocation
        'outputs': {'
            dnd': {'id': 'd5bb7278791be519',
                   'src': 'hda',
                   'uuid': '5e072cfd-45a5-4bf0-b17e-9e3ad7f72cc1'},
            'output': {'id': 'bf726321666e2d4e',
                       'src': 'hda',
                       'uuid': 'e4e15313-d355-4071-903e-24456c0c98cb'}}}

        '''
        print("..WRITING LOGS..")
        gi = self.galaxy_instance
        galaxy_pid = job_context['pid']
        state = job_context['state_dict']['state']
        # check that stdout is not empty
        n=0
        while n < 5:
            job = gi.jobs.show_job(galaxy_pid, full_details=True)
            print(f"Job: {job}")
            outputs = job['outputs']
            print(outputs)
            for _,dataset in outputs.items():
                dataset_info = gi.datasets.show_dataset(dataset['id'])
                print(f"Dataset info: {dataset_info}")
                provenance = gi.histories.show_dataset_provenance(history_id=get_history_id(gi,str(self.workspace)), dataset_id= dataset['id'])
                print(f"Provenance: {provenance}")
                if len(provenance['stdout']) != 0 or state == 'upload':
                    for std, file in dict(stderr=state + '_stderr', stdout=state + '_stdout').items():
                        print(f"writing galaxy {std}...in {file} for job {galaxy_pid}")
                        write_to_file(self.log_filenames_dict[file],
                                      'galaxy file name: ' + dataset_info['name'] + ' in job: ' + galaxy_pid + '\n' +
                                      provenance[std] + '\n' + dataset_info['name'] + '\n')
                    n = 5
                else:
                    print("waiting....")
                    time.sleep(1)


    def download_file(self,job_context):
        i = job_context['state_dict'].get('idx')
        result =  job_context['results'][i]
        remote_name = result.get('remote_name')
        download_path = os.path.join(self.LOCAL_WORKSPACE,self.workspace,result.get('file'))
        job_dir = os.path.join(self.LOCAL_WORKSPACE,self.workspace)
        file_ext = result.get('type')
        history_name = self.workspace
        self.connect()
        gi = self.galaxy_instance
        dataset = gi.histories.show_matching_datasets(history_id=get_history_id(gi,history_name), name_filter= remote_name)[0]
        download_url = dataset['download_url'] + '?to_ext=' + file_ext
        url = urljoin(gi.base_url, download_url)
        cmd = f"(nohup bash -c \"curl -o {download_path} {url} \" >/tmp/mtest </dev/null 2>/tmp/mtest.err & echo $!; wait $!; echo $? >> {job_dir}/$!.exit_status)"
        print(cmd)
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().rstrip()
        print(f"PID: {pid}")
        return pid

# GALAXY INIZIALIZATION

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


def convert_workflows_to_formly():
    workflow_files_list = os.listdir(os.path.join(ROOT,'biobarcoding/workflows'))

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
            w_id = workflow_id(gi, name)
            if w_id != 'there is no workflow named{}'.format(name):
                gi.workflows.delete_workflow(w_id)
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

# class getInputs(ToFormlyConverter):
#     def __init__(self):
#         super(getInputs, self).__init__(step_label)


class convertBooleanToolParameter(ToFormlyConverter):
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
    dictionary step_label: form_path
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

