from . import get_ngd_genus, push_taxonomy, RANKS
from ....io.taxa import TaxaTaskTools


def run():
	print("Synchronizing taxonomic hierarchy with Biota")

	print(">> Get data from ngd")
	ngd_genus = get_ngd_genus()

	print(">> Querying Biota")
	biota_df = TaxaTaskTools.biota_get_df()
	ranks = [_ for _ in RANKS[:-2] if _ in biota_df.columns]
	biota_df = biota_df[biota_df.genus.isin(ngd_genus)][ranks]

	biota_taxa = []

	def get_parent(lineage: list) -> str:
		_l = list(lineage[:-1])
		_l.reverse()
		return next((t for t in _l if t is not None), None)

	for i in range(len(biota_df.columns)):
		for _, r in biota_df[list(biota_df.columns[:i+1])].drop_duplicates().iterrows():
			biota_taxa.append({'canonicalName': r[-1], 'parent': get_parent(r), 'rank': r.keys()[-1]})

	push_taxonomy('biota', biota_taxa)
	return 'DONE'
