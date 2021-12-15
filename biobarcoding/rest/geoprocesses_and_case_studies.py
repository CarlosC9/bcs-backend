import json

from flask import g, request, Blueprint
from flask.views import MethodView
from sqlalchemy import and_

from .. import app_acronym
from ..authentication import n_session
from ..common import generate_json
from ..db_models.core import CaseStudy, CProcess, CProcessInstance, CProcessPort, Dataset, PortType, \
    CaseStudy2FunctionalObject
from . import make_simple_rest_crud, get_content, ResponseObject, Issue, IType, register_api, app_api_base
from ..services.geoprocesses import create_geoprocess_instance, submit_geoprocess_instances, create_geoprocess_instances


def load_case_study_item_to_delete(session, item_id):
    """
    Load a case study item to be deleted

    :param session:
    :param item_id:
    :return:
    """
    # Split the item_id
    item_id_parts = item_id.split(',')
    return session.query(CaseStudy2FunctionalObject).\
        filter(and_(CaseStudy2FunctionalObject.case_study_id == item_id_parts[0],
                    CaseStudy2FunctionalObject.functional_object_id == item_id_parts[1])).one()


bp_case_studies_fos, CaseStudiesFOSAPI = \
    make_simple_rest_crud(CaseStudy2FunctionalObject, "case_study_items",
                          alt_getters=dict(delete=load_case_study_item_to_delete))
bp_case_studies, CaseStudiesAPI = make_simple_rest_crud(CaseStudy, "case_studies")


# Geoprocesses
def custom_geoprocesses_filter(filter, session=None):
    """
    Obtain geoprocess in which datasets in "for_layers" can be used

    :param filter: A dictionary (str, dict)
    :return:
    """
    clauses = []
    if filter.get('for_layers'):
        ds_lst = ', '.join([f"{ds}" for ds in filter["for_layers"]["unary"]])
        sql_select = f"SELECT DISTINCT(g.id) " \
                     f"FROM {app_acronym}_processes g " \
                     f"JOIN {app_acronym}_processes_ports gp ON g.id=gp.process_id " \
                     f"WHERE gp.input AND " \
                     f"gp.port_type_id IN (" \
                     f"SELECT DISTINCT(port_type_id) FROM {app_acronym}_dataset_port_types d WHERE d.dataset_id in ({ds_lst})" \
                     f")"
        close_connection = session.transaction is None
        conn = session.connection()
        sql_result = conn.execute(sql_select)
        geoprocess_ids = [r[0] for r in sql_result]
        if len(geoprocess_ids) == 0:
            geoprocess_ids = [-1]
        clauses.append(CProcess.id.in_(geoprocess_ids))
        if close_connection:
            conn.close()
    return clauses


bp_geoprocesses, GeoprocessesAPI = make_simple_rest_crud(CProcess, "geoprocesses",
                                                         aux_filter=custom_geoprocesses_filter)
bp_geoprocess_port_types, GeoprocessPortTypesAPI = make_simple_rest_crud(PortType, "geoprocess_port_types")
bp_geoprocesses_ports, GeoprocessesPortsAPI = make_simple_rest_crud(CProcessPort, "geoprocesses_ports")


class GeoprocessInstancesAPI(MethodView):
    @n_session(read_only=True)
    def get(self, id_=None):
        """
        One geoprocess instance or a filtered & sorted list of them

        export API_BASE_URL=http://localhost:5000/api
        curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        curl --cookie app-cookies.txt "$API_BASE_URL/geo/processes/instances/"

        :param id_:
        :return:
        """
        def custom_instances_filter(filter, session=None):
            """
            case_studies

            sess.query(PortType).
            :param filter:
            :return:
            """
            clauses = []
            close_connection = session.transaction is None
            conn = session.connection()
            ids = set()
            avoid_no_entity = False

            if filter.get('for_geoprocesses'):  # May be used as input of a geoprocess (it could have been used already)
                _ = ', '.join(['\''+str(i)+'\'' for i in ids]) if ids else ''
                remaining_layers_condition = f" AND pi.id IN ({_})" if ids else ""
                lst = ', '.join([f"{ds}" for ds in filter["for_geoprocesses"]["unary"]])
                sql_select = f"SELECT DISTINCT(pi.id) " \
                             f"FROM {CProcessInstance.__tablename__} pi " \
                             f"JOIN {CProcess.__tablename__} g ON pi.instantiated_process_id=g.id " \
                             f"WHERE g.id IN ({lst})" \
                             f"{remaining_layers_condition}"
                sql_result = conn.execute(sql_select)
                _ = [r[0] for r in sql_result]
                if not _:
                    avoid_no_entity = True
                if ids:
                    ids.intersect(_)
                else:
                    ids = set(_)

            if filter.get("case_studies_"):
                _ = ', '.join(['\''+str(i)+'\'' for i in ids]) if ids else ''
                remaining_layers_condition = f" AND pi.id IN ({_})" if ids else ""
                lst = ', '.join([f"{ds}" for ds in filter["case_studies_"]["unary"]])
                sql_select = f"SELECT DISTINCT(pi.id) " \
                             f"FROM {CProcessInstance.__tablename__} pi " \
                             f"JOIN {CaseStudy2FunctionalObject.__tablename__} cs ON pi.id=cs.functional_object_id " \
                             f"WHERE cs.case_study_id IN ({lst}" \
                             f"{remaining_layers_condition}"
                sql_result = conn.execute(sql_select)
                ids = [r[0] for r in sql_result]
                _ = [r[0] for r in sql_result]
                if not _:
                    avoid_no_entity = True
                if ids:
                    ids.intersect(_)
                else:
                    ids = set(_)

            if avoid_no_entity and len(ids) == 0:
                ids = [-1]

            if ids and len(ids) > 0:
                clauses.append(CProcessInstance.id.in_(list(ids)))

            if close_connection:
                conn.close()

            return clauses

        def enhance_instance_json(instance: CProcessInstance):
            """
            Enhance the instance JSON with additional information

            :param instance:
            :return:
            """
            inputs = f""
            outputs = f""
            port_info = []
            for port in instance.ports:
                s = f"{port.port.name}: {port.dataset.name if port.dataset else '-'} ({port.dataset_id})"
                if port.port.input:
                    if inputs != "":
                        inputs += "; "
                    inputs += s

                else:
                    if outputs != "":
                        outputs += "; "
                    outputs += s
                port_info.append((port.dataset_id if port.dataset else None, port.port.input))
            _ = json.loads(generate_json(instance))
            for i in range(len(port_info)):
                _["ports"][i].update({"dataset_id": port_info[i][0], "input": port_info[i][1]})
            _["geoprocess"] = instance.instantiated_process.name
            _["geoprocess_id"] = instance.instantiated_process.id
            _["inputs"] = inputs
            _["outputs"] = outputs

            return _

        issues = []
        db_sess = g.n_session.db_session
        issues, content, count, status = get_content(db_sess, CProcessInstance, issues, id_,
                                                     aux_filter=custom_instances_filter)
        if isinstance(content, list):
            # Add "input_layers" summary and "output_layers" summary fields
            enh_content = []
            for inst in content:
                enh_content.append(enhance_instance_json(inst))
            content = enh_content
        else:
            content = enhance_instance_json(content)
        return ResponseObject(issues=issues, status=status, content=content, count=count).get_response()

    @n_session()
    def post(self):
        """
        Create a geoprocess instance

        Input JSON
        {
         "geoprocess": "<name (one of the names in geoprocesses.yaml)>"
         "inputs": {
            "<name1>": <layer_id1>,
            "<name2>": <layer_id2>
         },
        }

        export API_BASE_URL=http://localhost:5000/api
        curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"

        CREATE ONE INSTANCE:
        curl --cookie app-cookies.txt -i -XPOST "$API_BASE_URL/geo/processes/instances/" -H "Content-Type: application/json" -d @"/home/rnebot/GoogleDrive/AA_RIVUCAN/rvc-backend/tests/data_test/geoprocess_instance_post.json"
        curl --cookie app-cookies.txt -i -XPOST "http://localhost:5000/api/geoprocess/instances/" -d "{}"

        CREATE MANY INSTANCES:
        curl --cookie app-cookies.txt -i -XPOST "$API_BASE_URL/geo/processes/instances/" -H "Content-Type: application/json" -d @"/home/rnebot/GoogleDrive/AA_RIVUCAN/rvc-backend/tests/data_test/geoprocess_instances_update_post.json"
        curl --cookie app-cookies.txt -i -XPOST "$API_BASE_URL/geo/processes/instances/" -H "Content-Type: application/json" -d @"/home/rnebot/GoogleDrive/AA_RIVUCAN/rvc-backend/tests/data_test/geoprocess_instances_governance_post.json"

        :return:
        """
        session = g.n_session.db_session
        d = request.get_json()
        if "geoprocess" in d and "inputs" in d:
            i, issues = create_geoprocess_instance(session,
                                                   d["geoprocess"],
                                                   d["inputs"],
                                                   d.get("simulate", False),
                                                   d.get("case_studies", []),
                                                   d.get("schedule", False))
            if i is None:
                status = 400  # Bad request
            else:
                status = 200  # Ok
        elif "geoprocesses" in d and "datasets" in d:
            i, issues = create_geoprocess_instances(session,
                                                    d["geoprocesses"],
                                                    d["datasets"],
                                                    d.get("simulate", False),
                                                    d.get("case_studies", []),
                                                    d.get("schedule", False))
            if i is None:
                status = 400  # Bad request
            else:
                status = 200  # Ok
        else:
            status = 400  # Bad request
            issues = [Issue(IType.ERROR, f'Either "geoprocess" and "inputs" or "geoprocesses" and "datasets" '
                                         f'must be defined in the posted JSON.')]
            i = None
        return ResponseObject(issues=issues, status=status, content=i).get_response()

    @n_session()
    def put(self, id_=None):
        """
        Modify geoprocess instance. Two modifications:
        * State: from proposed to executing, from executing to executed correctly or incorrectly
        * Outputs: to specify the datasets for the output ports, once it has been executed

        export API_BASE_URL=http://localhost:5000/api
        curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        SCHEDULE a Geoprocess instance for execution (instance with ID 2 in this case)
        curl --cookie app-cookies.txt -i -XPUT "http://localhost:5000/api/geo/processes/instances/2" -H "Content-Type: application/json" -d '{"status": "scheduled"}'

        { "status": "...",
          "outputs": {
            "<name1>": <layer_id2>
          }
        }
        :param id_: ID of geoprocess instance
        :return:
        """
        issues = []
        status = 200  # By default, changed if there is something wrong
        content = None
        if id_:
            session = g.n_session.db_session
            d = request.get_json()
            i = session.query(CProcessInstance).get(id_)
            if not i or d is None:
                if not i:
                    issues.append(
                        Issue(IType.ERROR,
                              f"Geoprocess instance {id_} not found"))
                if not d:
                    issues.append(
                        Issue(IType.ERROR,
                              f"No JSON has been provided to PUT"))
                status = 400
            else:
                something_done = False
                if "status" in d:
                    req_status = d["status"]
                    # TODO TMP Remove "success" and "error" from the list below?
                    if i.status in ("candidate", "c", "cancelled", "success", "error") and req_status == "scheduled":
                        # Submit!!
                        job_ids = submit_geoprocess_instances(session, [i])
                        issues.append(
                            Issue(IType.INFO, f"Geoprocess instance {id_} is 'scheduled', with Job ID {job_ids[0]}"))
                        content = job_ids[0]
                        i.status = req_status
                        something_done = True
                    elif i.status == "scheduled" and req_status in ("cancelled", "success", "error"):
                        i.status = req_status
                        issues.append(Issue(IType.INFO, f"Geoprocess instance {id_} status is now {req_status}"))
                        something_done = True
                if "outputs" in d:
                    for name, layer_id in d["outputs"].items():
                        # Find port id
                        for port in i.ports:
                            if port.port.name.lower() == name.lower():
                                # Check if the dataset exists
                                ds = session.query(Dataset).get(layer_id)
                                if ds:
                                    port.dataset_id = layer_id  # Set dataset id
                                    issues.append(
                                        Issue(IType.INFO,
                                              f"Output port {name} of Geoprocess instance {id_} set to layer {layer_id}"))
                                else:
                                    issues.append(Issue(IType.ERROR, f"Geolayer with id {layer_id} not found"))
                                    status = 400
                                break
                if something_done:
                    session.commit()

        else:
            issues.append(Issue(IType.ERROR, 'Must specify ID of geoprocess instance to modify'))
            status = 400

        return ResponseObject(issues=issues, status=status, content=content).get_response()

    @n_session()
    def delete(self, id_=None):
        """
        Delete geoprocess instance (for statuses: "candidate", "cancelled" and "error")
        :param id_:
        :return:
        """
        issues = []
        status = 200  # By default, changed if there is something wrong
        content = None
        if id_:
            session = g.n_session.db_session
            i = session.query(CProcessInstance).get(id_)
            if not i:
                issues.append(
                    Issue(IType.INFO,
                          f"Geoprocess instance {id_} not found"))
                status = 400
            else:
                if i.status in ("candidate", "cancelled", "error"):
                    # TODO Cascade should work
                    # Remove process instance "ports", then the process instance
                    # session.query(PortInProcessInstance).filter(PortInProcessInstance.process_instance == i).delete()
                    session.delete(i)
                    issues.append(Issue(IType.INFO, f"Geoprocess instance {id_} deleted"))
                else:
                    issues.append(Issue(IType.ERROR, f"Geoprocess instance {id_} with status '{i.status} could not be deleted (only instances with statuses: 'candidate', 'cancelled', 'error')"))
                    status = 400
        else:
            issues.append(Issue(IType.ERROR, 'Must specify ID of geoprocess instance to delete'))
            status = 400

        return ResponseObject(issues=issues, status=status, content=content).get_response()


bp_geoprocess_instances = Blueprint(f'bp_geoprocess_instances', __name__)
entity_name = "geo/processes/instances"
_ = register_api(bp_geoprocess_instances, GeoprocessInstancesAPI, entity_name, f"{app_api_base}/{entity_name}/", "id_")
