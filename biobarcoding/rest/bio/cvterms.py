from flask import Blueprint
from .. import app_api_base
from ..crudie import CrudieAPI


def add_rules(bp=None) -> Blueprint:
    if not bp:
        from . import bp
    bp_view = CrudieAPI.as_view('api_cvterms', entity='cvterms')
    bp.add_url_rule(
        app_api_base + '/ontologies/terms',
        view_func=bp_view,
        methods=['GET']
    )
    bp.add_url_rule(
        app_api_base + '/ontologies/terms/',
        view_func=bp_view,
        methods=['GET']
    )
    bp.add_url_rule(
        app_api_base + '/ontologies/<int:cv_id>/terms',
        view_func=bp_view,
        methods=['GET']
    )
    bp.add_url_rule(
        app_api_base + '/ontologies/<int:cv_id>/terms/',
        view_func=bp_view,
        methods=['GET']
    )
    bp.add_url_rule(
        app_api_base + '/ontologies/terms/<int:cvterm_id>',
        view_func=bp_view,
        methods=['GET']
    )
    bp.add_url_rule(
        app_api_base + '/ontologies/<int:cv_id>/terms/<int:cvterm_id>',
        view_func=bp_view,
        methods=['GET']
    )
    return bp
