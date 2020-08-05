from bioblend import galaxy


def login(api_key, url):
    user_key = api_key
    gi = galaxy.GalaxyInstance(url=url, key=user_key)
    return gi


def library_list(gi):
    libraries = gi.libraries.get_libraries()
    return libraries


def library_id(gi,libname: 'library name to search'):
    libraries = gi.libraries.get_libraries()
    lib = next(item for item in libraries if item["name"] == libname)
    return lib['id']


def workflow_list(gi):
    workflows = gi.workflows.get_workflows()
    return workflows


def workflow_id(gi, name: 'library name to search'):
    workflows = gi.workflows.get_workflows()
    work = next(item for item in workflows if item["name"] == name)
    return work['id']


def dataset_list(gi):
    datasets = gi.datasets.get_datasets()
    return datasets


def dataset_id(gi, name: 'libraty name to search'):
    datasets = dataset_list(gi)
    dat = next(item for item in datasets if item["name"] == name)
    return dat['id']


def create_library(gi, library_name):
    '''
    look for a library or creates a new library if it donst exists.
    Returns the library id of the dirst library finded with that name
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



def create_input(step, source,id):
    '''
    Create a rorkflow input
    :param step:  step input
    :param source: 'hdda' (history); 'lbda' (library)
    :param id: dataset id
    :return: Input dictionary
    '''
    datamap = dict()
    datamap[step] = {'id': id, 'src': source}
    return datamap


def run_workflow(gi , name: 'str', input_path: 'str',input_name, *, step_index: 'str' ='0', history_name: 'str' ='Test_History'):
    '''
    Directly run a workflow from a local file
    #1 Create New History
    #2 Upload Dataset
    #3 create the new input usinf upload_file tool in the History just created
    #4 Invoke the Workflow in the new History -> this Generates: id (como WorkflowInvocation), workflow_id, history_id, y cada step de ese Workflow invocation tendrá un Job ID
    '''
    h_id = gi.histories.create_history(name=history_name)['id']
    d_id = gi.tools.upload_file(input_path, h_id, filename= input_name)['outputs'][0]['id']
    w_id = workflow_id(gi, name)
    dataset = {'src': 'hda', 'id': d_id}
    invocation = gi.workflows.invoke_workflow(w_id,
                                              inputs={step_index: dataset},
                                              history_id=h_id,
                                              inputs_by='step_index',
                                              )
    return invocation

def list_invocation_results(gi,invocation_id):
    '''
    :param gi: Galaxy Instance
    :param invocation_id: Workgflow Invocation Id
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




def download_result(gi, results, path):
    '''

    :param gi: Galaxy Instance
    :param invocation_id:Workgflow Invocation Id
    :param name: list of dictionaties of datasets
    :return: download results
    '''
    if isinstance(results,list):
        for r in results:
            gi.datasets.download_dataset(r['id'], file_path=path)
    else:
        print(results)



def show_jobs():
    pass








