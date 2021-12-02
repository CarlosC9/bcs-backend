import json

from flask import g, request, Blueprint
from flask.views import MethodView

from . import ResponseObject, get_content, Issue, IType, register_api, app_api_base, \
    make_simple_rest_crud
from ..authentication import n_session
from ..db_models.hierarchies import Hierarchy, HierarchyNode
from ..services.hierarchies import get_hierarchy, create_or_update_hierarchy, \
    create_update_or_delete_hierarchy_nodes, get_just_hierarchy

bp_hierarchy_nodes, HierarchyNodesAPI = make_simple_rest_crud(HierarchyNode, "hierarchy_nodes")


class HierarchiesAPI(MethodView):
    @n_session(read_only=True)
    def get(self, id_=None):
        """
        Hierarchies

        export API_BASE_URL=http://localhost:5000/api
        curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        curl --cookie app-cookies.txt "$API_BASE_URL/hierarchies/"

        :param id_:
        :return:
        """
        issues = []
        db_sess = g.n_session.db_session
        if id_ is None:
            issues, content, count, status = get_content(db_sess, Hierarchy, issues, id_=id_)
        else:
            content = get_hierarchy(db_sess, id_, "flat")
            if content is None:
                status = 400
            else:
                status = 200
                # content = json.dumps(content)
        return ResponseObject(issues=issues, status=status, content=content, count=count).get_response()

    @n_session()
    def post(self):
        """
        Create a hierarchy

        Input JSON, based on the dict:
          dict([uuid], [type_id], name, attributes, nodes=[], levels=[])

        where levels:
          dict([uuid], name)

        where nodes:
          dict([uuid], name, level=<l>, parent=<n>, reference_node=<n>, attributes=dict())

        export API_BASE_URL=http://localhost:5000/api
        curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        curl --cookie app-cookies.txt -i -XPOST "$API_BASE_URL/hierarchies/" -H "Content-Type: application/json" -d @"/home/rnebot/GoogleDrive/AA_RIVUCAN/rvc-backend/tests/data_test/hierarchy_post.json"
        curl --cookie app-cookies.txt -i -XPOST "$API_BASE_URL/hierarchies/" -H "Content-Type: application/json" -d "{}"

        :return:
        """
        session = g.n_session.db_session
        d = request.get_json()
        h = create_or_update_hierarchy(session,
                                       d.get("uuid"), d.get("name", "undefined"), d.get("type_id", 0),
                                       d.get("attributes", {}))
        if h is None:
            issue = Issue(IType.ERROR, 'Could not create Hierarchy')
            status = 400  # Bad request
        else:
            if d.get("nodes"):
                # Create nodes
                create_update_or_delete_hierarchy_nodes(session, h, d["nodes"])
            session.commit()
            issue = Issue(IType.INFO, f'Hierarchy created correctly. (ID: {h.id})')
            status = 200  # Ok

        return ResponseObject(issues=[issue], status=status, content=h).get_response()

    @n_session()
    def put(self, id_=None):
        """
        Modify hierarchy

        export API_BASE_URL=http://localhost:5000/api
        curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"
        curl --cookie app-cookies.txt -i -XPUT "$API_BASE_URL/hierarchies/2" -H "Content-Type: application/json" -d

        :param id_: ID of hierarchy
        :return:
        """
        issues = []
        status = 200  # By default, changed if there is something wrong
        content = None
        if id_:
            session = g.n_session.db_session
            d = request.get_json()
            h = get_just_hierarchy(session, id_)
            if not h or d is None:
                if not h:
                    issues.append(Issue(IType.ERROR, f"Hierarchy {id_} not found"))
                if d is None:
                    issues.append(Issue(IType.ERROR, f"No JSON has been provided to PUT"))
                status = 400
            else:
                h = create_or_update_hierarchy(session,
                                               d.get("uuid"), d.get("name", "undefined"), d.get("type_id", 0),
                                               d.get("attributes", {}))
                if d.get("nodes"):
                    # Create nodes
                    create_update_or_delete_hierarchy_nodes(session, h, d["nodes"])

                session.commit()
                content = h
        else:
            issues.append(Issue(IType.ERROR, 'Must specify ID of hierarchy to modify'))
            status = 400

        return ResponseObject(issues=issues, status=status, content=content).get_response()

    @n_session()
    def delete(self, id_=None):
        """

        :param id_:
        :return:
        """
        issues = []
        status = 200  # By default, changed if there is something wrong
        content = None
        return ResponseObject(issues=issues, status=status, content=content).get_response()


bp_hierarchies = Blueprint(f'bp_hierarchies', __name__)
entity_name = "hierarchies"
_ = register_api(bp_hierarchies, HierarchiesAPI, entity_name, f"{app_api_base}/{entity_name}/", "id_")

