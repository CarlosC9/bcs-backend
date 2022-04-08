from flask import Blueprint
from biobarcoding.rest.crudie import CrudieAPI

bp = Blueprint('bp_bio', __name__)

from .bos import add_rules
_bp = add_rules(bp)
from .metadata import add_rules
_bp = add_rules(_bp)
from .cvterms import add_rules
_bp = add_rules(_bp)

bp_bio = _bp