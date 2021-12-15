"""
REST interface to manage JOBS API
"""
import json

from dotted.collection import DottedDict
from flask import Blueprint
from flask import request, make_response, Response, jsonify, g
from flask.views import MethodView
from sqlalchemy import and_

from ..authentication import n_session
from ..common import generate_json
from ..common.helpers import is_integer
from ..db_models import DBSession
from ..db_models.jobs import ProcessInComputeResource
from ..jobs import JobManagementAPI
from ..jobs.process_adaptor import ProcessAdaptorFactory
from . import app_api_base, register_api, Job, ComputeResource, Process, ResponseObject, \
    SocketService, decode_request_params, IType, Issue

bp_jobs = Blueprint('jobs', __name__)


def get_job_executability(job, session):
    # Check if execution can proceed, or it must be contained
    max_running_jobs = job.resource.jm_params.get("max_running_jobs", 0)
    if max_running_jobs > 0:
        states = ','.join([f"'{s}'" for s in ["preparing_workspace",
                                              "export",
                                              "transfer_data_to_resource",
                                              "submit",
                                              "wait_until_execution_starts",
                                              "wait_for_execution_end",
                                              "transfer_data_from_resource",
                                              "store_result_in_backend",
                                              "cleanup"]])
        # First, find the count of running jobs in the resource
        sql_count = f"select count(*) as cnt from jobs_jobs j where j.resource_id = {job.resource_id} and status in ({states})"
        close_connection = session.transaction is None
        conn = session.connection()
        sql_result = conn.execute(sql_count)
        cont = sql_result.first()
        cont = cont[0]
        if cont < max_running_jobs:
            # Check if the job has priority (job id)
            s = f"SELECT count(*) as CNT " \
                f"FROM " \
                f"(select id from jobs_jobs j where j.resource_id = {job.resource_id} and status = 'created' " \
                f" order by id limit {max_running_jobs - cont}) as a " \
                f"WHERE id = {job.id}"
            sql_result = conn.execute(s)
            res = sql_result.first()
            res = res[0]
            result = "resource_ok_must_wait" if res == 0 else "resource_ok_can_execute"
        else:
            result = "resource_ko"
        if close_connection:
            conn.close()
    else:
        result = "resource_ok_can_execute"
    return result


# Jobs REST API
class JobAPI(MethodView):
    """
    Job Resource
    """
    page: int = None
    page_size: int = None
    decorators = []  # Add decorators: identity, function execution permissions, logging, etc.

    @n_session(read_only=True)
    def get(self, job_id=None):
        # return "<h1 style='color:blue'>Hello JOBS!</h1>"
        session = g.n_session.db_session
        r = ResponseObject()
        final_statuses = ['success', 'cancelled', 'error']

        status = request.args.get("status")
        if job_id is None:
            query = session.query(Job)
            if status:
                if status == "done":
                    query = query.filter(Job.status.in_(final_statuses))
                elif status == "todo":
                    query = query.filter(Job.status.notin_(final_statuses))
            self.__check_data(request.args)
            if self.page and self.page_size:
                query = query.offset((self.page - 1) * self.page_size).limit(self.page_size)
            r.content = query.all()
        else:
            # Return detail and status of single job
            query = session.query(Job).filter(Job.id == job_id)
            if status:
                if status == "done":
                    query = query.filter(Job.status.in_(final_statuses))
                elif status == "todo":
                    query = query.filter(Job.status.notin_(final_statuses))
            job = query.first()
            if job:
                r.content = job
                if job.status == "created":
                    _ = json.loads(generate_json(job))
                    _["executability"] = get_job_executability(job, session)
                    r.content = _
            else:
                r.status = 400
                r.issues = [Issue(IType.ERROR, f"Job with ID {job_id} not found")]
                r.content = None

        return r.get_response()

    @n_session(read_only=True)
    def post(self):
        """
        curl -i -XPOST http://localhost:5000/api/jobs/ --data-urlencode "{}"
        curl -i -XPOST http://localhost:5000/api/jobs/ -H "Content-Type: application/json" -d @"/home/paula/Documentos/App/backend/tests/data_test/test_job_req.json"
        :return:
        """
        # Submit new Job
        msg = f'POST {request.path}\nPosting job'

        # Get Identity ID
        identity_id = g.n_session.identity_id
        if identity_id is None:
            return Response("User not authorized", status=401)

        # Start session
        session = g.n_session.db_session

        # Start JSON for processing
        d = DottedDict()
        uncoded_req = request.get_json()
        req = decode_request_params(uncoded_req.get('params'))
        # Load resource and process
        print(req)
        in_dict = DottedDict(req)

        if is_integer(in_dict.resource_id):
            resource = session.query(ComputeResource).get(in_dict.resource_id)
        else:
            resource = session.query(ComputeResource).filter(ComputeResource.uuid == in_dict.resource_id).first()
        if is_integer(in_dict.process_id):
            process = session.query(Process).get(in_dict.process_id)
        else:
            process = session.query(Process).filter(Process.uuid == in_dict.process_id).first()
        process_in_resource = session.query(ProcessInComputeResource).filter(
            and_(ProcessInComputeResource.process_id == process.id,
                 ProcessInComputeResource.resource_id == resource.id)).first()
        process_params = in_dict.process_params
        """
        Input JSON
        {
            "resource_id": ...
            "process_id": ...
            "process_params": {
            }
            "credentials": {
            }
        """
        # Prepare "job_context" for Celery tasks
        # "job_id": ...,
        # "endpoint_url": ...
        # "resource": {
        #     "name": "",
        #     "jm_type": "",
        #     "jm_location": {
        #     },
        #     "jm_credentials": {
        #     }
        # }
        # "process": {
        #     "name": "",
        #     "inputs": {
        #     }
        # }

        # PROCESS
        d.process = DottedDict()
        d.status = "created"
        d.process.inputs = process_params
        d.process.name = process_in_resource.native_process_id
        # RESOURCE
        d.resource = DottedDict()
        d.resource.name = resource.name
        d.identity_id = identity_id
        d.resource.jm_type = resource.jm_type.name
        d.resource.jm_location = resource.jm_location
        d.resource.jm_credentials = resource.jm_credentials if "credentials" not in in_dict else in_dict.credentials
        process_adaptor = ProcessAdaptorFactory().get(d.resource.jm_type, in_dict.process_id)
        d = process_adaptor.adapt_job_context(d)

        outputs = [r.to_json() for r in d.results]
        # Create Job database object
        job = Job()
        job.resource = resource
        job.process = process
        job.status = d.status
        job.identity_id = identity_id
        job.inputs = json.dumps(process_params.to_json())
        job.outputs = json.dumps(outputs)
        session.add(job)
        session.commit()
        d.job_id = job.id
        DBSession.remove()
        # Submit job to Celery

        JobManagementAPI().submit(d.to_json())
        # Return
        job_dict = job.__dict__
        response_object = {
            'status': 'success',
            'message': f"Job with ID {job.id} queued",
            'job': dict((k, job_dict[k]) for k in job_dict if k not in ["_sa_instance_state", "uuid"])
        }

        return make_response(jsonify(response_object)), 200

    @n_session()
    def delete(self, job_id):
        # Delete Job. Only statuses: created, finished unsuccessfully, cancelled
        db = g.n_session.db_session
        job = db.query(Job).filter(Job.id == job_id).first()
        if job:
            # TODO TMP Remove "cancelling" !!! (use only for tests)
            can_hard_delete = job.status in ("cancelled", "error", "cancelling")
            can_soft_delete = job.status == "success"
            status = 200
            issue_type = IType.INFO
            if can_hard_delete:
                # Check if there is a CProcessInstance referring to it. Then, put that reference to None before delete
                from ..db_models.core import CProcessInstance
                inst = db.query(CProcessInstance).filter(CProcessInstance.job_id == job_id).first()
                if inst:
                    inst.job_id = None
                    if inst.status == "scheduled":
                        inst.status = "candidate"  # Return to candidate, because Job execution was incomplete
                db.delete(job)
                msg = f'DELETE {request.path}\nHard-deleting job {job_id} with current status: {job.status}'
            elif can_soft_delete:
                if not job.deleted:
                    job.deleted = True
                    msg = f'DELETE {request.path}\nSoft-deleting job {job_id} with current status: {job.status}'
                else:
                    status = 400
                    msg = f'DELETE {request.path}\nJob {job_id} with status: {job.status}, already Soft-deleted'
                    issue_type = IType.ERROR
            else:
                msg = f'DELETE {request.path}\nJob {job_id} with status: {job.status}. Cancel first with PUT status=cancelling'
                issue_type = IType.ERROR
                status = 400
        else:
            status = 400
            msg = f'DELETE {request.path}\nJob {job_id} does not exist'
            issue_type = IType.ERROR

        r = ResponseObject(issues=[Issue(issue_type, msg)],
                           status=status)
        return r.get_response()

    @n_session()
    def put(self, job_id):
        # Update job? For now, only allow changing status
        # (also, return "status", because cannot change "cancelled" status)
        db = g.n_session.db_session
        req = request.get_json()
        job = db.query(Job).filter(Job.id == job_id).first()
        if job:
            status = 200
            i_type = IType.INFO
            job_status = job.status
            if job.status != req['status']:
                if job.status not in ("cancelling", "cancelled", "success", "error") or \
                        (job.status == "cancelling" and req["status"] == "cancelled"):
                    msg = f'PUT {request.path}\nModifying job {job_id} status from "{job.status}" to "{req["status"]}"'
                    job.status = req['status']
                    job_status = job.status
                    db.add(job)
                    if SocketService.instance:
                        job_dict = job.__dict__
                        socket_job_dict = dict(
                            (k, job_dict[k]) for k in job_dict if k not in ["_sa_instance_state", "uuid"])
                        SocketService.instance.emit_process_status(job.identity_id, socket_job_dict)
                else:
                    msg = f'PUT {request.path}\nRequest modify job {job_id} ' \
                          f'status to "{req["status"]}" but it is "{job.status}"'
            else:
                msg = f'PUT {request.path}\nRequested status for job {job_id}, "{job.status}", is already set'
        else:
            i_type = IType.ERROR
            status = 400
            msg = f'PUT {request.path}\nRequested status for job {job_id}, but Job does not exist'
            job_status = None
        print(msg)
        r = ResponseObject(issues=[Issue(i_type, msg)],
                           status=status,
                           content=dict(status=job_status, requested_status=req["status"]))
        return r.get_response()

    def __check_data(self, data):
        self.page = int(data.get(self.page, 1))
        self.page_size = int(data.get(self.page_size, 1000000))


register_api(bp_jobs, JobAPI, "jobs", f"{app_api_base}/jobs/", pk="job_id")
