import hashlib
import traceback
from io import BytesIO

import canonicaljson
import itertools
import json
import subprocess

import yaml
from dotted.collection import DottedDict
from flask import g
from sqlalchemy import func
from typing import List, Tuple

from .. import get_global_configuration_variable
from ..common.helpers import is_integer
from ..db_models import DBSession
from ..db_models.core import CProcessInstance, CProcess, PortInProcessInstance, CProcessPort, Dataset, PortType, \
    CaseStudy, CaseStudy2FunctionalObject
from ..db_models.jobs import ComputeResource, JobManagementType, Process, Job
from ..jobs import JobManagementAPI
from ..jobs.ssh_process_adaptors import SSHProcessAdaptor
from ..rest import Issue, IType, logger
from ..services.files import put_file_contents


def get_geoprocess_definitions(resource):
    port = 22 if "port" not in resource.jm_location else resource.jm_location["port"]
    s = f"/usr/bin/ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -p {port} {resource.jm_credentials['username']}@{resource.jm_location['host']} \"geoprocess definitions\""
    p = subprocess.Popen(s, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, err = p.communicate()
    rc = p.returncode
    return rc, output, err


def submit_geoprocess_instances(session, geoprocess_instances):
    def expand_params(params, process: CProcess) -> DottedDict:
        _ = DottedDict()
        _.geoprocess = params.geoprocess
        _.inputs = []
        for name, layer_id in params.inputs.items():
            # Find port in "process.inputs"
            for port in process.ports:
                if port.input and port.name.lower() == name.lower():
                    extension = port.attributes.get("format", "gpkg").lower()
                    break
            _.inputs.append(DottedDict(dict(name=name, layer_id=layer_id, type=extension, file=f"{name}.{extension}")))
        _.outputs = []
        for name, layer_id in params.outputs.items():
            # Find port in "process.inputs"
            for port in process.ports:
                if not port.input and port.name.lower() == name.lower():
                    extension = port.attributes.get("format", "gpkg").lower()
                    if extension == "gpkg":
                        content_type = "application/geopackage+sqlite3"
                    else:
                        content_type = "application/binary"
                    break
            tmp = port.attributes.copy()
            tmp.update(dict(port_id=port.id))
            _.outputs.append(DottedDict(dict(name=name, type=extension, content_type=content_type, attributes=tmp,
                                             file=f"{name}.{extension}")))

        return _

    def submit(ident, resource: ComputeResource, instance: CProcessInstance):
        """
        Submit "geoprocess instance" to Jobs

        :param ident:
        :param resource:
        :param instance:
        :return:
        """
        # Load persistent objects
        process = session.query(Process).filter(Process.name == "geoprocess").first()  # Load process "geoprocess"
        params = DottedDict(instance.params)
        process_def = session.query(CProcess).filter(
            func.lower(CProcess.name) == params.geoprocess.lower().strip()).first()
        # Generate and upload params.yml to FilesAPI
        params = expand_params(params, process_def)
        params_file_path = f"geoprocess/instance_{instance.id}/params.yml"
        put_file_contents(session,
                          f"{params_file_path}.content",
                          BytesIO(yaml.dump(params.to_python(), encoding='ascii')),
                          "text/yaml")

        # Create "job context" without process adaptor, because it will be executed by an SSH job manager
        d = DottedDict()
        d.status = "created"
        d.endpoint_url = get_global_configuration_variable("ENDPOINT_URL")
        # Resource
        d.resource = DottedDict()
        d.resource.name = resource.name
        d.resource.jm_type = resource.jm_type.name
        d.resource.jm_location = resource.jm_location
        d.resource.jm_credentials = resource.jm_credentials
        # Process
        d.process = DottedDict()
        d.process.name = f"geoprocess {params.geoprocess}"
        d.process.inputs = DottedDict()
        d.process.inputs.parameters = DottedDict()
        d.process.inputs.parameters[SSHProcessAdaptor.SCRIPT_KEY] = ["geoprocess"]
        d.process.inputs.parameters[SSHProcessAdaptor.SCRIPT_FILES_KEY] = []
        d.process.inputs.parameters[SSHProcessAdaptor.SCRIPT_PARAMS_KEY] = "execute ."
        d.process.callbacks = DottedDict()
        d.process.callbacks.endpoint = f"/geo/processes/instances/{instance.id}"
        d.process.callbacks.success_status = "success"
        d.process.callbacks.error_status = "error"
        d.process.callbacks.cancelled_status = "cancelled"
        d.results = []
        # params.yml "->"
        d.process.inputs.data = []
        d.process.inputs.data.append(dict(remote_name="params.yml",
                                          type=None,
                                          object_type={"files": params_file_path},
                                          selection=0))
        # Input layers "->"
        for i_ in params.inputs:
            d.process.inputs.data.append(dict(remote_name=i_.file,
                                              type=i_.type,
                                              object_type={"geo": "layers"},
                                              selection=i_.layer_id))

        # Output layers "<-"
        for o in params.outputs:
            _ = o.attributes.to_python()
            # Pass information necessary to associate the resulting dataset with the instance port
            _.update(dict(instance_id=instance.id))
            d.results.append(dict(remote_name=o.file,
                                  file=o.file,
                                  type=o.type,
                                  content_type=o.content_type,
                                  autoimport=True,
                                  name=f"capa automática",
                                  metadata=dict(attributes=_),
                                  object_type={"geo": "layers"}))

        # Create Job database object
        job = Job()
        job.resource = resource
        job.process = process
        job.status = d.status
        job.identity_id = ident
        job.inputs = d.process.inputs.to_python()
        job.outputs = d.results.to_python()
        session.add(job)
        session.flush()
        instance.job_id = job.id
        session.commit()
        d.job_id = job.id

        # # Generate CURL for test purposes
        # cookies_file_path = get_global_configuration_variable("COOKIES_FILE_PATH")
        # file_dict = d.results[0]
        # file_dict["metadata"]["job_id"] = job.id
        # endpoint = get_global_configuration_variable("ENDPOINT_URL")
        # curl_cmd = f'curl -s --cookie-jar {cookies_file_path} --cookie {cookies_file_path} ' \
        #            f'-F \'layer_file=@{file_dict["file"]};type={file_dict["content_type"]}\' ' \
        #            f'-F \'metadata={json.dumps(file_dict["metadata"])};type=application/json\' ' \
        #            f'\'{endpoint}{app_api_base}/geo/layers/\''
        # print(curl_cmd)

        # Submit job to Celery
        JobManagementAPI().submit(d.to_json())

        return job.id

    # ------------------------------------------------------------------------------------------------------------------

    identity_id = g.n_session.identity_id
    if identity_id is None:
        return None

    # ComputeResources, "geoprocess" enabled
    ssh_resources = session.query(ComputeResource).join(ComputeResource.jm_type).\
        filter(JobManagementType.name == "ssh").all()
    resources = []
    for r in ssh_resources:
        rc, output, err = get_geoprocess_definitions(r)
        if rc == 0:
            resources.append(r)

    # Submit to the different resources
    job_ids = []
    for i, inst in enumerate(geoprocess_instances):
        job_id = submit(identity_id, resources[i % len(resources)], inst)
        if job_id:
            inst.status = "scheduled"
        job_ids.append(job_id)

    return job_ids


def create_geoprocess_instance(session, geoprocess_name, ports, simulate=False, case_studies=None, schedule=False) -> \
        Tuple[CProcessInstance, List[Issue]]:
    """
    Create a persistent geoprocess instance

    :param session: SQLAlchemy session
    :param geoprocess_name: Geoprocess name
    :param ports: List of what goes into each of the input ports
    :param simulate: True: do not persist
    :param case_studies: To which case studies associate the instance and the instance results
    :param schedule: True: schedule for execution
    :return:
    """
    if is_integer(geoprocess_name):
        process_def = session.query(CProcess).get(int(geoprocess_name))
    elif isinstance(geoprocess_name, str):
        process_def = session.query(CProcess).filter(func.lower(CProcess.name) == geoprocess_name.lower().strip()).first()
    elif isinstance(geoprocess_name, CProcess):
        process_def = geoprocess_name
    else:
        process_def = None
    issues = []
    if process_def:
        inst_desc = f"Proc: '{process_def.name}' [{process_def.id}]"
        params = dict(geoprocess=process_def.name, inputs={}, outputs={})
        i = CProcessInstance()
        i.instantiated_process = process_def
        i.status = "candidate"
        rollback = False
        # Ports
        outputs = []
        first_input = True
        for port in process_def.ports:  # For each port in the Geoprocess Definition
            if not port.input:  # Outputs, just placeholder
                p = PortInProcessInstance()
                p.process_instance = i
                p.port = port
                session.add(p)
                params["outputs"][port.name] = None
                outputs.append(p)
            else:  # Inputs, fully specified
                layer_id = ports.get(port.name, ports.get(port, 0))
                if layer_id is not None:
                    if isinstance(layer_id, Dataset):
                        dataset = layer_id
                        layer_id = dataset.id
                    else:
                        dataset = session.query(Dataset).get(layer_id)
                        if dataset is None:
                            rollback = True
                            issues.append(Issue(IType.ERROR, f'Could not find Dataset {layer_id} for port {port.name}'))

                    if first_input:
                        inst_desc += "("
                        first_input = False
                    else:
                        inst_desc += f", "
                    inst_desc += f"{port.name}: {dataset.name+' ['+str(dataset.id)+']' if dataset else ''}"
                    p = PortInProcessInstance()
                    p.process_instance = i
                    p.port = port
                    p.dataset = dataset
                    session.add(p)
                    params["inputs"][port.name] = layer_id
                else:
                    # ERROR - Could not find layer id for input port name {port.name}
                    rollback = True
                    issues.append(Issue(IType.ERROR, f'Could not find specification of Dataset for port {port.name}'))
        inst_desc += ")"
        # Unique instance?
        _ = canonicaljson.encode_canonical_json(dict(process=i.instantiated_process.id, params=params))
        _ = hashlib.sha256(_).hexdigest()
        equiv_instances = session.query(CProcessInstance).filter(CProcessInstance.hash_of_canonical == _).all()
        if len(equiv_instances) > 0:
            # ERROR - Instance already exists!
            issues.append(Issue(IType.ERROR, f'Another instance of geoprocess "{process_def.name}" with exactly the'
                                             f' same input Datasets was created previously '
                                             f'(with ID: {equiv_instances[0].id})'))
            rollback = True
        else:
            i.hash_of_canonical = _

        # Case studies
        _ = set(case_studies)
        cs_objects = session.query(CaseStudy).filter(CaseStudy.id.in_(_)).all()
        if len(cs_objects) != len(_):
            # ERROR - Did not find all specified case studies
            issues.append(Issue(IType.ERROR, f'Could not find some of the case studies '
                                             f'specified to associate the instance with'))
            rollback = True
        else:
            # TODO Pass case studies to each output port (they would inherit the membership automatically)
            pass

        i.params = params

        # Final processing
        if rollback:
            session.rollback()
            return None, issues
        else:
            session.add(i)
            for cs in cs_objects:
                _ = CaseStudy2FunctionalObject()
                _.case_study = cs
                _.functional_object = i
            if not simulate:
                session.commit()
                if schedule:
                    submit_geoprocess_instances(session, [i])
            else:
                session.rollback()

            return i, [Issue(IType.INFO, f'Geoprocess instance created '
                                         f'{"(simulated)" if simulate else ""} '
                                         f'correctly: {inst_desc} (ID: {i.id})')]
    else:
        # ERROR - Geoprocess not defined
        return None, [Issue(IType.ERROR, f'Geoprocess to instantiate not defined correctly')]


def create_geoprocess_instances(session, geoprocess_names, datasets_to_consider, simulate=False, case_studies=None, schedule=False):
    """
    Creation of geoprocess instances looking for matches

    :param session:
    :param geoprocess_names: List of geoprocess names
    :param datasets_to_consider:
    :param simulate: True: do not create instances. Just return the
    :param case_studies: To which case studies associate the instance and the instance results
    :param schedule: True: schedule for execution
    :return:
    """
    geoprocesses = {}
    result = []
    issues = []
    for ds in datasets_to_consider:
        if is_integer(ds):
            ds = session.query(Dataset).get(int(ds))
        types = [t.port_type_id for t in ds.types]  # Port types
        if len(types) == 0:
            issues.append(Issue(IType.WARNING, f'The dataset {ds.id} does not have associated Types'))
            continue
        for geoprocess_name in geoprocess_names:
            if is_integer(geoprocess_name):
                gp = session.query(CProcess).get(int(geoprocess_name))
                geoprocesses[gp.name.lower()] = gp
                name = gp.name
            elif isinstance(geoprocess_name, str):
                name = geoprocess_name
                if name.lower() not in geoprocesses:
                    geoprocesses[name.lower()] = session.query(CProcess).filter(
                        func.lower(CProcess.name) == name.lower().strip()).first()
            elif isinstance(geoprocess_name, CProcess):
                name = geoprocess_name.name
                geoprocesses[name] = geoprocess_name
            for port in geoprocesses[name.lower()].ports:
                if port.input and port.port_type_id in types:
                    # Match!
                    any_zero = False
                    other_ports = []
                    for port2 in geoprocesses[name.lower()].ports:
                        if port2 != port and port2.input:
                            # Matching datasets
                            datasets = session.query(Dataset).join(Dataset.types).filter(
                                PortType.id == port2.port_type_id).all()
                            if len(datasets) == 0:
                                any_zero = True
                                break
                            else:
                                other_ports.append([(port2, ds2) for ds2 in datasets])

                    if not any_zero:
                        # Generate geoprocess instances
                        already_created_combinations = set()
                        for i, combination in enumerate(list(itertools.product(*other_ports))):
                            t_ = tuple(sorted([f"{t[0].id}:{t[1].id}" for t in combination] + [f"{port.id}:{ds.id}"]))
                            if t_ in already_created_combinations:
                                logger.debug(f'Combination {t_} already created')
                                continue
                            already_created_combinations.add(t_)
                            ports = {t[0]: t[1] for t in combination}
                            ports[port] = ds
                            try:
                                res, _ = create_geoprocess_instance(session, geoprocesses[name.lower()], ports,
                                                                    simulate=simulate, case_studies=case_studies,
                                                                    schedule=False)
                            except Exception as e:
                                traceback.print_exc()
                                res = None
                                _ = [Issue(IType.ERROR, f'Exception {e} creating geoprocess instance')]

                            issues.extend(_)
                            if res:
                                result.append(res)
                        # Schedule
                        if not simulate and schedule:
                            submit_geoprocess_instances(session, result)
    return result, issues


def update_geoprocesses():
    # For each "compute resource with ssh", execute: ssh -p 8022 root@localhost "geoprocess definitions"
    #  If a JSON is returned, synchronize with CProcess definitions
    resources = DBSession.query(ComputeResource).join(ComputeResource.jm_type).\
        filter(JobManagementType.name == "ssh").all()
    for r in resources:
        rc, output, err = get_geoprocess_definitions(r)
        if rc == 0:
            # Read Geoprocess Definitions
            geoprocesses = json.loads(output)
            session = DBSession()
            # Initialize PortTypes
            port_types = {}  # str -> PortType
            for geoprocess in geoprocesses:
                for j, ports_list in enumerate([geoprocess["inputs"], geoprocess["outputs"]]):
                    for json_port in ports_list:
                        pt = json_port.get("port_type", None)
                        if pt:
                            if pt.lower() not in port_types:
                                port_type = session.query(PortType).filter(
                                    func.lower(PortType.name) == pt.lower().strip()).first()
                                if not port_type:
                                    port_type = PortType()
                                    port_type.name = pt
                                    port_type.attributes = {}
                                    session.add(port_type)
                                port_types[pt.lower()] = port_type

            # Geoprocesses
            for geoprocess in geoprocesses:
                i = session.query(CProcess).filter(func.lower(CProcess.name) == geoprocess["name"].lower().strip()).first()
                if not i:
                    i = CProcess()
                    # Set attributes
                    i.name = geoprocess["name"]
                    i.process_type = geoprocess["type"]
                    i.attributes = {}
                    session.add(i)
                    ports = []
                else:
                    ports = session.query(CProcessPort).filter(CProcessPort.process_id == i.id).all()
                # Update attributes directly
                excludes = {"name", "type", "inputs", "outputs"}
                i.attributes.update({k: geoprocess[k] for k in set(geoprocess.keys()).difference(excludes)})

                # Port names
                p_names = set([p.name.lower() for p in ports])
                # Now, update geoprocess ports
                for j, ports_list in enumerate([geoprocess["inputs"], geoprocess["outputs"]]):
                    for json_port in ports_list:
                        if json_port["name"].lower() not in p_names:
                            port = CProcessPort()
                            port.name = json_port["name"]
                            port.input = j == 0
                            port.process = i
                            port.attributes = {}
                            session.add(port)
                        else:
                            for port in ports:
                                if port.name.lower() == json_port["name"].lower():
                                    break
                        # Update port Attributes
                        excludes = {"name"}
                        port.attributes.update({k: json_port[k] for k in set(json_port.keys()).difference(excludes)})
                        # Update port Type
                        pt = json_port.get("port_type", None)
                        if pt:
                            port_type = port_types[pt.lower()]
                        else:
                            port_type = None
                        port.port_type = port_type

            session.commit()
        else:
            print(f"Error: {err}")
