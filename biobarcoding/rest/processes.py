from biobarcoding.rest import make_simple_rest_crud
from biobarcoding.db_models.jobs import Process

bp_processes, ProcessesAPI = make_simple_rest_crud(Process, "processes")