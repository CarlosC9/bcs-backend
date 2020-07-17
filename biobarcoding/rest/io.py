from flask import Blueprint

bp_io = Blueprint('io', __name__)

from biobarcoding.services import FileAPI

file = FileAPI.as_view('file_api')

bp_io.add_url_rule(
    '/io/export',
    view_func=file,
    methods=['GET']
)
