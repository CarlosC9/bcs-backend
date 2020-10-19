"""
REST interface to manage JOBS API
"""
from dotted.collection import DottedDict, DottedCollection
from flask import Blueprint
from flask import request, make_response, jsonify
from flask.views import MethodView
import json
from biobarcoding.db_models import DBSession
from biobarcoding.jobs import JobManagementAPI
from biobarcoding.rest import bcs_api_base, register_api, Job, ComputeResource, Process

bp_jobs = Blueprint('jobs', __name__)


# Jobs REST API

class JobAPI(MethodView):
    """
    Job Resource
    """

    decorators = []  # Add decorators: identity, function execution permissions, logging, etc.

    def get(self, job_id):
        return "<h1 style='color:blue'>Hello JOBS!</h1>"
        if job_id is None:
            tmp = JobManagementAPI.list(filter)
            return make_response(jsonify(tmp)), 200
        else:
            # Return detail and status of single job
            tmp = JobManagementAPI.get(job_id)
            return make_response(jsonify(tmp)), 200

    def post(self):
        """
        curl -i -XPOST http://localhost:5000/api/jobs/ --data-urlencode "{}"
        :return:
        """
        # Submit new Job
        msg = f'POST {request.path}\nPosting job'
        # Start session
        session = DBSession()
        # Start JSON for processing
        d = DottedDict()
        d.endpoint_url = ""
        # Load resource and process
        in_dict = DottedCollection.load_json("{}")
        # resource = session.query(ComputeResource).get(in_dict.resource.id)
        process = session.query(Process).filter(Process.name == "pd-1.0").first()
        d.resource = DottedDict()
        d.resource.name = ""  # resource.name
        # Create Job database object
        job = Job()
        job.status = "created"
        job.resource = None
        job.process = process
        session.add(job)
        session.commit()
        d.job_id = job.id
        DBSession.remove()

        # Prepare JSON for Celery
        # "resource": {
        #     "name": "",
        #     "proc_arch": "",
        #     "operating_system": "",
        #     "jm_type": "",
        #     "jm_location": {
        #     },
        #     "jm_credentials": {
        #
        #     }
        # }
        # "workflow": {
        #     "name": "",
        #     "id": "",
        #     "inputs": {
        #     }
        # }

        # Submit job to Celery
        JobManagementAPI().submit(d.to_json())

        # Return
        response_object = {
            'status': 'success',
            'message': f"Job with ID {job.id} queued"
        }

        return make_response(jsonify(response_object)), 200

    def delete(self, job_id):
        # Cancel Job
        msg = f'DELETE {request.path}\nDeleting job {id}'

        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject)), 200

    def put(self, job_id):
        # Update job? What would be the utility
        msg = f'PUT {request.path}\nModifying job {id}'
        print(msg)
        self._check_data()

        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject)), 200


register_api(bp_jobs, JobAPI, "jobs", f"{bcs_api_base}/jobs/", "job_id")
