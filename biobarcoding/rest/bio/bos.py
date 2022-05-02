from flask import Blueprint
from .. import app_api_base
from ..crudie import CrudieAPI


def add_rules(bp=None) -> Blueprint:
    if not bp:
        from . import bp
    for bos in ['sequences', 'alignments', 'phylotrees', 'discriminant_matrices', 'blasts', 'supermatrices']:
        bp_view = CrudieAPI.as_view('api_'+bos, entity=bos)
        bp.add_url_rule(
            app_api_base + f'/bos/{bos}',
            view_func=bp_view,
            methods=['GET', 'POST', 'DELETE']
        )
        bp.add_url_rule(
            app_api_base + f'/bos/{bos}/',
            view_func=bp_view,
            methods=['GET', 'POST', 'DELETE']
        )
        bp.add_url_rule(
            app_api_base + f'/bos/{bos}/<string:id>',
            view_func=bp_view,
            methods=['GET', 'PUT', 'DELETE']
        )
        bp.add_url_rule(
            app_api_base + f'/bos/{bos}.<string:format>',
            view_func=bp_view,
            methods=['GET']
        )
        bp.add_url_rule(
            app_api_base + f'/bos/{bos}/.<string:format>',
            view_func=bp_view,
            methods=['GET']
        )
        bp.add_url_rule(
            app_api_base + f'/bos/{bos}/<string:id>.<string:format>',
            view_func=bp_view,
            methods=['GET']
        )
    return bp