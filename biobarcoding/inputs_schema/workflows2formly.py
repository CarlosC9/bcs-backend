import json
import os

from ..common import ROOT


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


class convertBooleanToolParameter(ToFormlyConverter):
    def __init__(self, step_label):
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
    def __init__(self, step_label):
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
    def __init__(self, step_label):
        super(converterTextToolParameter, self).__init__(step_label)

    def convert(self, input):
        form = self.conversion(input)
        form['type'] = 'textarea'
        return form


class converterConditional(ToFormlyConverter):
    def __init__(self, step_label):
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


class inputCreation:

    def __init__(self, workflow_path):
        self.workflow = workflow_path

    @property
    def workflow(self):
        return self.__workflow

    @workflow.setter
    def workflow(self, workflow_path):
        with open(workflow_path) as jsonfile:
            self.__workflow = json.load(jsonfile)

    def __get_inputs(self):
        input = []
        workflow = self.workflow
        steps = workflow.get("steps")
        for _, step in steps.items():
            if len(step.get("input_connections")) == 0:
                for step_input in step["inputs"]:
                    input.append({"name": step_input.get("name"),
                                  "bo_type": step_input.get("bo_type", "no specified"),
                                  # TODO DONDE ESTA ESTA INFORMACIÓN
                                  "bo_format": step_input.get("bo_format", "no specified")
                                  }
                                 )

        return input

    def convert(self):
        inputs = self.__get_inputs()
        inputs_fieldGroups = list()
        for input in inputs:
            form = dict()
            form["key"] = "input"
            form["templateOptions"] = {"label": "input"}
            form['fieldGroup'] = [{
                "key": "remote_name",
                "type": "input",
                "templateOptions": {
                    "label": "input name"
                },
                "defaultValue": input["name"]
            },
                {
                    "key": "type",
                    "type": "input",
                    "templateOptions": {
                        "label": "type"
                    },
                    "defaultValue": input["bo_format"]
                },
                {
                    "key": "bo_type",
                    "type": "input",
                    "templateOptions": {
                        "label": "bo type"
                    },
                    "defaultValue": input["bo_type"]
                }]
            inputs_fieldGroups.append(form)
        # a esto le tengo que hacer el append de los parámetros de los algoritmos
        return inputs_fieldGroups


def convert_workflows_to_formly():
    # workflow_files_list = os.listdir(os.path.join(ROOT,'biobarcoding/workflows'))

    wfdict1 = {'steps': {'clustalw': ROOT + '/biobarcoding/inputs_schema/clustalw_galaxy.json',
                         'phyml': ROOT + '/biobarcoding/inputs_schema/phyml_galaxy.json'}
        ,
               'fname': 'clustalw_phyml_formly.json',
               'workflow_path': '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/biobarcoding/workflows/Galaxy-Workflow-ClustalW-PhyMl.ga'
               }
    wfdict2 = {'steps':
                   [{'clustalw': ROOT + '/biobarcoding/inputs_schema/clustalw_galaxy.json'}
                    ],
               'fname': 'clustalw_formly.json',
               'workflow_path': '/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/biobarcoding/workflows/Galaxy-Workflow-MSA_ClustalW.ga'
               }

    path = ROOT + '/biobarcoding/inputs_schema/'
    lwdict = [wfdict1, wfdict2]
    for wfdict in lwdict:
        if wfdict['fname'] not in os.listdir(path):
            wfdict['fname'] = path + wfdict['fname']
            convertToFormly(wfdict['steps'], wfdict['workflow_path'], wfdict['fname'])


def convertToFormly(wf_steps, path_to_workflow, newpath):
    '''
    dictionary step_label: form_path
    '''
    creator = inputCreation(path_to_workflow)
    fieldGroup = creator.convert()
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
    file = open(newpath, 'w')
    file.write(formly_json)
    file.close()
