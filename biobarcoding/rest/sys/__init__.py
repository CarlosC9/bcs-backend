from flask import Blueprint

bp = Blueprint('bp_sys', __name__)

from .sysadmin import add_rules
_bp = add_rules(bp)
from .formats import add_rules
_bp = add_rules(bp)

bp_sys = _bp