from flask import Blueprint

bp_msa = Blueprint('msa', __name__)

from biobarcoding.services import AlignAPI

msa = AlignAPI.as_view('msa_api')
bp_msa.add_url_rule(
    '/bo/msa/',
    view_func=msa,
    methods=['GET','POST']
)
bp_msa.add_url_rule(
    '/bo/msa/<int:msa_id>',
    view_func=msa,
    methods=['GET','PUT','DELETE']
)

from biobarcoding.services import AlignFeatAPI

msa_feat = AlignFeatAPI.as_view('msa_feat_api')
bp_msa.add_url_rule(
    '/bo/msa/<msa_id>/feature/',
    view_func=msa_feat,
    methods=['GET','POST']
)
bp_msa.add_url_rule(
    '/bo/msa/<int:msa_id>/feature/<int:cmt_id>',
    view_func=msa_feat,
    methods=['GET','PUT','DELETE']
)
