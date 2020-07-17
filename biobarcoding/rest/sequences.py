from flask import Blueprint

bp_seq = Blueprint('seq', __name__)

from biobarcoding.services import SeqAPI

seq = SeqAPI.as_view('seq_api')
bp_seq.add_url_rule(
    '/bo/sequence',
    view_func=seq,
    methods=['GET','POST']
)
bp_seq.add_url_rule(
    '/bo/sequence/<int:id>',
    view_func=seq,
    methods=['GET','PUT','DELETE']
)

from biobarcoding.services import SeqFeatAPI

seq_feat = SeqFeatAPI.as_view('seq_feat_api')
bp_seq.add_url_rule(
    '/bo/sequence/<seq_id>/feature/',
    view_func=seq_feat,
    methods=['GET','POST']
)
bp_seq.add_url_rule(
    '/bo/sequence/<int:seq_id>/feature/<int:cmt_id>',
    view_func=seq_feat,
    methods=['GET','PUT','DELETE']
)
