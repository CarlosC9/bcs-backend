from flask import Blueprint

bp_jobqueue = Blueprint('jobqueue', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

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
bp_jobqueue.add_url_rule(
    '/job/queue/<int:job_type>/<int:job_id>',
    view_func=job_queue,
    methods=['GET','PUT','DELETE']
)
bp_jobqueue.add_url_rule(
    '/job/queue/<int:job_type>',
    view_func=job_queue,
    methods=['GET']
)
bp_jobqueue.add_url_rule(
    '/job/queue/',
    view_func=job_queue,
    methods=['GET']
)
