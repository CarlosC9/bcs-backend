from biobarcoding.rest import make_simple_rest_crud,ResponseObject
from biobarcoding.db_models.jobs import Process, ComputeResource, ProcessInComputeResource
from flask.views import MethodView
from flask import request, make_response, jsonify, send_file

bp_processes, ProcessesAPI = make_simple_rest_crud(Process, "processes")
bp_resources, ResourcesAPI = make_simple_rest_crud(ComputeResource, "resources")