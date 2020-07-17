from flask import Blueprint

bp_meta = Blueprint('metadata', __name__)

from biobarcoding.services import OntologyAPI

onto = OntologyAPI.as_view('onto_api')
bp_meta.add_url_rule(
    '/bo/ontology',
    view_func=onto,
    methods=['GET','POST']
)
bp_meta.add_url_rule(
    '/bo/ontology/<int:id>',
    view_func=onto,
    methods=['GET','DELETE']
)

from biobarcoding.services import TaxonomyAPI

taxon = TaxonomyAPI.as_view('taxon_api')
bp_meta.add_url_rule(
    '/bo/organism',
    view_func=taxon,
    methods=['GET','POST']
)
bp_meta.add_url_rule(
    '/bo/organism/<int:id>',
    view_func=taxon,
    methods=['GET','DELETE']
)

from biobarcoding.services import AnalysisAPI

ansis = AnalysisAPI.as_view('ansis_api')
bp_meta.add_url_rule(
    '/bo/analysis',
    view_func=ansis,
    methods=['GET','POST']
)
bp_meta.add_url_rule(
    '/bo/analysis/<int:id>',
    view_func=ansis,
    methods=['GET','DELETE']
)
