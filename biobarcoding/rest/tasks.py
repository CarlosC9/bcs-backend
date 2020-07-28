"""
REST for management of Celery Tasks
"""

from flask import Blueprint

bp_tasks = Blueprint('tasks', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


# Task API

class TaskAPI(MethodView):
    """
    Task Resource
    """

    def post(self, type):
        msg = f'POST {request.path}\nCreating task {args}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def put(self, id):
        msg = f'PUT {request.path}\nModifying task {id}'
        print(msg)
        self._check_data()

        responseObject = {
        'status': 'success',
        'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, id):
        msg = f'DELETE {request.path}\nDeleting task {id}'
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


task = TaskAPI.as_view('task_api')
bp_tasks.add_url_rule(
    bcs_api_base + '/task/run/<task_type>',
    view_func=task,
    methods=['POST']
)
bp_tasks.add_url_rule(
    bcs_api_base + '/task/run/<int:task_id>',
    view_func=task,
    methods=['PUT','DELETE']
)


# Task Queue API

class TaskQueueAPI(MethodView):
    """
    TaskQueue Resource
    """
    def get(self, type=None, id=None):
        msg = f'GET {request.path}\nGetting task queue {type}/{id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def put(self, type, id):
        msg = f'PUT {request.path}\nCreating task queue {type}/{id}'
        print(msg)
        self._check_data()

        responseObject = {
            'status': 'success',
            'message': msg
        }
        return make_response(jsonify(responseObject)), 200


    def delete(self, type, id):
        msg = f'DELETE {request.path}\nDeleting task queue {type}/{id}'
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


task_queue = TaskQueueAPI.as_view('task_queue_api')
bp_tasks.add_url_rule(
    bcs_api_base + '/task/queue/<int:task_type>/<int:task_id>',
    view_func=task_queue,
    methods=['GET','PUT','DELETE']
)
bp_tasks.add_url_rule(
    bcs_api_base + '/task/queue/<int:task_type>',
    view_func=task_queue,
    methods=['GET']
)
bp_tasks.add_url_rule(
    bcs_api_base + '/task/queue/',
    view_func=task_queue,
    methods=['GET']
)
