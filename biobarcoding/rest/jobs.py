"""
REST interface to manage JOBS API
"""
from dotted.collection import DottedDict, DottedList, DottedJSONEncoder, DottedCollection
from flask import Blueprint
from flask import request, make_response, jsonify, Response, g
from flask.views import MethodView
import json

from sqlalchemy import and_
from alchemyjsonschema import SchemaFactory, StructuralWalker

from biobarcoding.common.helpers import is_integer
from biobarcoding.db_models import DBSession
from biobarcoding.db_models.jobs import ProcessInComputeResource
from biobarcoding.jobs import JobManagementAPI
from biobarcoding.jobs.process_adaptor import ProcessAdaptorFactory
from biobarcoding.rest import bcs_api_base, register_api, Job, ComputeResource, Process, ResponseObject
from biobarcoding.authentication import bcs_session
from biobarcoding.rest import make_simple_rest_crud

bp_jobs = Blueprint('jobs', __name__)


# Jobs REST API
class JobAPI(MethodView):
    """
    Job Resource
    """
    page: int = None
    page_size: int = None
    decorators = []  # Add decorators: identity, function execution permissions, logging, etc.

    @bcs_session(read_only=True)
    def get(self, job_id = None):
        # return "<h1 style='color:blue'>Hello JOBS!</h1>"
        db = g.bcs_session.db_session
        r = ResponseObject()
        if job_id is None:
            query = db.query(Job)
            self.__check_data(request.args)
            if self.page and self.page_size:
                query = query.offset((self.page - 1) * self.page_size).limit(self.page_size)
            r.content = query.all()
        else:
            # Return detail and status of single job
            r.content = db.query(Job).filter(Job.id == job_id).first()
        return r.get_response()

    def post(self):
        """
        curl -i -XPOST http://localhost:5000/api/jobs/ --data-urlencode "{}"
        curl -i -XPOST http://localhost:5000/api/jobs/ -H "Content-Type: application/json" -d @"/home/paula/Documentos/NEXTGENDEM/bcs/bcs-backend/tests/data_test/test_job_req.json"
        :return:
        """
        # Submit new Job
        msg = f'POST {request.path}\nPosting job'
        # Start session
        session = DBSession()
        # Start JSON for processing
        d = DottedDict()
        req = request.get_json()
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
        process_in_resource = session.query(ProcessInComputeResource).filter(and_(ProcessInComputeResource.process_id==process.id, ProcessInComputeResource.resource_id==resource.id)).first()
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
        d.resource.jm_type = resource.jm_type.name
        d.resource.jm_location = resource.jm_location
        d.resource.jm_credentials = resource.jm_credentials if "credentials" not in in_dict else in_dict.credentials

        # Create Job database object
        job = Job()
        job.resource = resource
        job.process = process
        job.status = d.status
        job.inputs = json.dumps(process_params.to_json())
        session.add(job)
        session.commit()
        d.job_id = job.id
        DBSession.remove()

        process_adaptor = ProcessAdaptorFactory().get(d.resource.jm_type, in_dict.process_id)
        d = process_adaptor.adapt_job_context(d)

        # Submit job to Celery
        JobManagementAPI().submit(d.to_json())

        # Return
        response_object = {
            'status': 'success',
            'message': f"Job with ID {job.id} queued"
        }

        return make_response(jsonify(response_object)), 200

    @bcs_session()
    def delete(self, job_id):
        # Cancel Job
        msg = f'DELETE {request.path}\nDeleting job {id}'

        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject)), 200

    @bcs_session()
    def put(self, job_id):
        # Update job? What would be the utility
        msg = f'PUT {request.path}\nModifying job {job_id}'
        print(msg)
        db = g.bcs_session.db_session
        # r = ResponseObject()
        req = request.get_json()
        job = db.query(Job).filter(Job.id == job_id).first()
        job.status = req['status']
        db.add(job)
        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject)), 200

    def __check_data(self, data):
        self.page = int(data.get(self.page, 1))
        self.page_size = int(data.get(self.page_size, 1000000))


register_api(bp_jobs, JobAPI, "jobs", f"{bcs_api_base}/jobs/", pk ="job_id")
