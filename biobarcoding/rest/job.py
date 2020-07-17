from flask import Blueprint

bp_job = Blueprint('job', __name__)

from biobarcoding.services import JobAPI

job = JobAPI.as_view('job_api')
bp_job.add_url_rule(
    '/job/run/<job_type>',
    view_func=job,
    methods=['POST']
)
bp_job.add_url_rule(
    '/job/run/<int:job_id>',
    view_func=job,
    methods=['PUT','DELETE']
)

from biobarcoding.services import JobQueueAPI

job_queue = JobQueueAPI.as_view('job_queue_api')
bp_job.add_url_rule(
    '/job/queue/<int:job_type>/<int:job_id>',
    view_func=job_queue,
    methods=['GET','PUT','DELETE']
)
bp_job.add_url_rule(
    '/job/queue/<int:job_type>',
    view_func=job_queue,
    methods=['GET']
)
bp_job.add_url_rule(
    '/job/queue/',
    view_func=job_queue,
    methods=['GET']
)
