from flask import Blueprint

bp_phylo = Blueprint('phylo', __name__)

from biobarcoding.services import PhyloAPI

phylo = PhyloAPI.as_view('phylo_api')
bp_phylo.add_url_rule(
    '/bo/phylo/',
    view_func=phylo,
    methods=['GET','POST']
)
bp_phylo.add_url_rule(
    '/bo/phylo/<int:phylo_id>',
    view_func=phylo,
    methods=['GET','PUT','DELETE']
)

from biobarcoding.services import PhyloFeatAPI

phylo_feat = PhyloFeatAPI.as_view('phylo_feat_api')
bp_phylo.add_url_rule(
    '/bo/phylo/<phylo_id>/feature/',
    view_func=phylo_feat,
    methods=['GET','POST']
)
bp_phylo.add_url_rule(
    '/bo/phylo/<int:phylo_id>/feature/<int:cmt_id>',
    view_func=phylo,
    methods=['GET','PUT','DELETE']
)
