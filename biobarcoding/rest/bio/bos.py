from flask import Blueprint
from .. import app_api_base
from ..crudie import CrudieAPI


# TODO: add discriminant-matrices, blasts, supermatrices, collections
def add_rules(bp=None) -> Blueprint:
    if not bp:
        from . import bp
    for bos in ['sequences', 'alignments', 'phylotrees']:
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
        # bp.add_url_rule(
        #     app_api_base + f'/bos/{bos}/<int:seq_id>/features/',
        #     view_func=bp_view,
        #     methods=['GET', 'POST']
        # )
        # bp.add_url_rule(
        #     app_api_base + f'/bos/{bos}/<int:seq_id>/features/<int:cmt_id>',
        #     view_func=bp_view,
        #     methods=['GET', 'PUT', 'DELETE']
        # )
    return bp