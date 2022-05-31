from flask import Blueprint
from .. import app_api_base
from ..crudie import CrudieAPI


def add_rules(bp=None) -> Blueprint:
    if not bp:
        from . import bp
    for entity in ['status_checkers']:
        bp_view = CrudieAPI.as_view('api_'+entity, entity=entity)
        bp.add_url_rule(
            app_api_base + f'/sys/{entity}',
            view_func=bp_view,
            methods=['GET', 'POST', 'PUT', 'DELETE']
        )
        bp.add_url_rule(
            app_api_base + f'/sys/{entity}/',
            view_func=bp_view,
            methods=['GET', 'POST', 'PUT', 'DELETE']
        )
        bp.add_url_rule(
            app_api_base + f'/sys/{entity}/<string:id>',
            view_func=bp_view,
            methods=['GET', 'PUT', 'DELETE']
        )
        bp.add_url_rule(
            app_api_base + f'/sys/{entity}.<string:format>',
            view_func=bp_view,
            methods=['GET']
        )
        bp.add_url_rule(
            app_api_base + f'/sys/{entity}/.<string:format>',
            view_func=bp_view,
            methods=['GET']
        )
        bp.add_url_rule(
            app_api_base + f'/sys/{entity}/<string:id>.<string:format>',
            view_func=bp_view,
            methods=['GET']
        )
    return bp
