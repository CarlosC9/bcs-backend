from operator import and_

from biobarcoding.authentication import bcs_session
from biobarcoding.rest import make_simple_rest_crud, ResponseObject, register_api
from biobarcoding.db_models.jobs import Process, ComputeResource, ProcessInComputeResource
from flask.views import MethodView
from flask import Response, Blueprint, g, request
from flask import request, make_response, jsonify, send_file

bp_processes, ProcessesAPI = make_simple_rest_crud(Process, "processes")
# bp_resources, ResourcesAPI = make_simple_rest_crud(ComputeResource, "resources")
# bp_resourcesinprocess,ResourcesInProcessAPI = make_simple_rest_crud(ProcessInComputeResource,"process/{id}/resource")

bcs_api_base = "/api"  # Base for all RESTful calls
bcs_gui_base = "/gui"  # Base for the Angular2 GUI
bcs_external_gui_base = "/gui_external"  # Base for the Angular2 GUI when loaded from URL


# TODO  no funciona
class ResourcesInProcessAPI(MethodView):
    page: int = None
    page_size: int = None

    @bcs_session(read_only=True)
    def get(self, pid, _id=None):
        db = g.bcs_session.db_session
        r = ResponseObject()
        if _id is None:
            # TODO corregir
            query = db.query(ComputeResource).join(ProcessInComputeResource).filter(Process.id == pid)
            self.__check_data(request.args)
            if self.page and self.page_size:
                query = query.offset((self.page - 1) * self.page_size).limit(self.page_size)
            r.content = query.all()
        else:
            # Detail
            r.content = db.query(ComputeResource).join(ProcessInComputeResource).filter(
                Process.id == pid).filter(ComputeResource.id == _id).first()

        return r.get_response()

    @bcs_session(authr=None)
    def post(self, pid):  # Create
        db = g.bcs_session.db_session
        r = ResponseObject()
        t = request.json
        # TODO check process id.. etc
        # s = ProcessInComputeResource.Schema().load(t, instance=ProcessInComputeResource())
        # db.add(s)
        return r.get_response()

    @bcs_session(authr=None)
    def put(self, pid, _id):  # Update (total or partial)
        db = g.bcs_session.db_session
        r = ResponseObject()
        t = request.json
        # s = db.query(ProcessInComputeResource).filter(
        #     and_(ProcessInComputeResource.process_id == pid, ProcessInComputeResource.resource_id == _id)).first()
        # s = ProcessInComputeResource.Schema().load(t, instance=s)
        # db.add(s)
        return r.get_response()

    @bcs_session(authr=None)
    def delete(self, pid, _id):  # Delete
        db = g.bcs_session.db_session
        r = ResponseObject()
        s = db.query(ProcessInComputeResource).filter(
            and_(ProcessInComputeResource.process_id == pid, ProcessInComputeResource.resource_id == _id)).first()
        db.delete(s)
        return r.get_response()

    def __check_data(self, data):
        self.page = int(data.get(self.page, 1))
        self.page_size = int(data.get(self.page_size, 1000000))


resource_process_view = ResourcesInProcessAPI.as_view('api_processes_resources')
bp_processes.add_url_rule(
    f"{bcs_api_base}/processes/<int:pid>/resources/<int:_id>/",
    view_func=resource_process_view,
    methods=['GET']
)
bp_processes.add_url_rule(
    f"{bcs_api_base}/processes/<int:pid>/resources/",
    view_func=resource_process_view,
    methods=['GET']
)
