from flask import Blueprint
from .. import app_api_base, ResponseObject
from ...services.main import get_service


def add_rules(bp=None) -> Blueprint:
    if not bp:
        from . import bp

    @bp.route(app_api_base + f'/sys/formats/<string:entity>')
    def get(entity):
        return ResponseObject(content=get_service(entity).formats).get_response()

    return bp
