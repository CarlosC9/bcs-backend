"""
REST interface to manage JOBS API
"""

from flask import Blueprint
from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base, register_api

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
            # Return ALL jobs (how to receive filter?) visible by the user
            pass
        else:
            # Return status of single job
            pass

    def post(self):
        # Submit new Job
        msg = f'POST {request.path}\nPosting job'
        responseObject = {
        'status': 'success',
        'message': msg
        }

        return make_response(jsonify(responseObject)), 200

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
