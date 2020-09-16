from bioblend import galaxy


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


def set_parameters(gi,workflow_id,nstep, param_name, new_value):
    workflow_info = gi.workflows.show_workflow(workflow_id)
    step = workflow_info['steps'][nstep]
    params = dict()
    for i in range(len(workflow_info['steps'])):
        params[str(i)] = workflow_info['steps'][str(i)]['tool_inputs']
    params[nstep][param_name] = new_value
    # Some Parse:
    for i in params:
        params[i] = {k: v.strip('"') for k, v in params[i] .items()}
        params[i].pop('input', None) # check if this is ok for every step (maybe not good for step = '0') or only necessary for workflows
    return params


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
    # TODO delete history after sucesfull download
    # gi.histories.delete_history(history_id)










