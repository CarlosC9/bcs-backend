from flask import Blueprint
from .. import app_api_base
from ..crudie import CrudieAPI


def add_rules(bp=None) -> Blueprint:
    if not bp:
        from . import bp
    for metadata in ['taxonomies', 'organisms', 'ontologies', 'analyses', 'individuals', 'collections']:
        bp_view = CrudieAPI.as_view('api_'+metadata, entity=metadata)
        bp.add_url_rule(
            app_api_base + f'/{metadata}',
            view_func=bp_view,
            methods=['GET', 'POST', 'DELETE']
        )
        bp.add_url_rule(
            app_api_base + f'/{metadata}/',
            view_func=bp_view,
            methods=['GET', 'POST', 'DELETE']
        )
        bp.add_url_rule(
            app_api_base + f'/{metadata}/<string:id>',
            view_func=bp_view,
            methods=['GET', 'PUT', 'DELETE']
        )
        bp.add_url_rule(
            app_api_base + f'/{metadata}.<string:format>',
            view_func=bp_view,
            methods=['GET']
        )
        bp.add_url_rule(
            app_api_base + f'/{metadata}/.<string:format>',
            view_func=bp_view,
            methods=['GET']
        )
        bp.add_url_rule(
            app_api_base + f'/{metadata}/<string:id>.<string:format>',
            view_func=bp_view,
            methods=['GET']
        )
    return bp
