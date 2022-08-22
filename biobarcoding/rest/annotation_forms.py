from flask import Blueprint
from . import app_api_base
from .crudie import CrudieAPI

bp_annotations = Blueprint('bp_annotations', __name__)

for item_type in ['form_templates', 'form_fields', 'form_relationships']:
    bp_view = CrudieAPI.as_view('api_annotation_'+item_type, entity=item_type)
    bp_annotations.add_url_rule(
        app_api_base + f'/annotation_{item_type}/',
        view_func=bp_view,
        methods=['GET', 'POST']
    )
    bp_annotations.add_url_rule(
        app_api_base + f'/annotation_{item_type}/<int:id>',
        view_func=bp_view,
        methods=['GET', 'PUT', 'DELETE']
    )
    bp_annotations.add_url_rule(
        app_api_base + f'/annotation_{item_type}/cvterm:<int:cvterm_id>',
        view_func=bp_view,
        methods=['GET', 'PUT', 'DELETE']
    )
    bp_annotations.add_url_rule(
        app_api_base + f'/annotation_{item_type}/<string:db>:<string:dbxref>',
        view_func=bp_view,
        methods=['GET', 'PUT', 'DELETE']
    )

bp_view = CrudieAPI.as_view('api_annotations', entity='annotations')
bp_annotations.add_url_rule(
    app_api_base + '/annotations',
    view_func=bp_view,
    methods=['GET', 'POST', 'DELETE']
)
bp_annotations.add_url_rule(
    app_api_base + '/annotations/',
    view_func=bp_view,
    methods=['GET', 'POST', 'DELETE']
)
bp_annotations.add_url_rule(
    app_api_base + '/annotations/<int:id>',
    view_func=bp_view,
    methods=['GET', 'PUT', 'DELETE']
)
bp_annotations.add_url_rule(
    app_api_base + '/annotations/<string:object_uuid>',
    view_func=bp_view,
    methods=['GET', 'POST', 'PUT', 'DELETE']
)

bp_view = CrudieAPI.as_view('api_relationships', entity='relationships')
bp_annotations.add_url_rule(
    app_api_base + '/relationships',
    view_func=bp_view,
    methods=['GET', 'POST', 'DELETE']
)
bp_annotations.add_url_rule(
    app_api_base + '/relationships/',
    view_func=bp_view,
    methods=['GET', 'POST', 'DELETE']
)
bp_annotations.add_url_rule(
    app_api_base + '/relationships/<int:id>',
    view_func=bp_view,
    methods=['GET', 'PUT', 'DELETE']
)
