"""
REST interface to manage JOBS API
"""

from flask import Blueprint

bp_jobs = Blueprint('jobs', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


# Job API

class JobAPI(MethodView):
    """
    Job Resource
    """

    def post(self, type):
        msg = f'POST {request.path}\nCreating job {args}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def put(self, id):
        msg = f'PUT {request.path}\nModifying job {id}'
        print(msg)
        self._check_data()

        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting job {id}'
        print(msg)
        self._check_data()

        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def _check_data(self):

        post_data = request.get_json()
        print(f'JSON data: {post_data}')


job = JobAPI.as_view('job_api')
bp_jobs.add_url_rule(
    bcs_api_base + '/job/run/<job_type>',
    view_func=job,
    methods=['POST']
)
bp_jobs.add_url_rule(
    bcs_api_base + '/job/run/<int:job_id>',
    view_func=job,
    methods=['PUT','DELETE']
)


# Job Queue API

class JobQueueAPI(MethodView):
    """
    JobQueue Resource
    """
    def get(self, type=None, id=None):
        msg = f'GET {request.path}\nGetting job queue {type}/{id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def put(self, type, id):
        msg = f'PUT {request.path}\nCreating job queue {type}/{id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, type, id):
        msg = f'DELETE {request.path}\nDeleting job queue {type}/{id}'
        print(msg)
        self._check_data()

        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def _check_data(self):

        post_data = request.get_json()
        print(f'JSON data: {post_data}')


job_queue = JobQueueAPI.as_view('job_queue_api')
bp_jobs.add_url_rule(
    bcs_api_base + '/job/queue/<int:job_type>/<int:job_id>',
    view_func=job_queue,
    methods=['GET','PUT','DELETE']
)
bp_jobs.add_url_rule(
    bcs_api_base + '/job/queue/<int:job_type>',
    view_func=job_queue,
    methods=['GET']
)
bp_jobs.add_url_rule(
    bcs_api_base + '/job/queue/',
    view_func=job_queue,
    methods=['GET']
)
