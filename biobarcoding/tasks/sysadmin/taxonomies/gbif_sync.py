from threading import Thread
from pygbif import species

from . import get_ngd_genus, push_taxonomy


def run(**kwargs):
	print("Synchronizing taxonomic hierarchy with GBIF")

	print(">> Get data from ngd")
	ngd_genus = get_ngd_genus()

	def get_gbif_key(org) -> int:
		return org.get('key', org.get('usageKey'))

	print(">> Querying GBIF")
	gbif_taxa, gbif_missing = [], []
	_th = Thread()
	for genus in ngd_genus:
		g_org = species.name_backbone(genus, kingdom='plantae', rank='genus')
		g_key = get_gbif_key(g_org)
		if not g_key:
			print('no key for', g_org)
			gbif_missing.append(g_org)
			continue
		g_parents = species.name_usage(g_key, data='parents')
		if g_parents:
			g_org['parent'] = g_parents[-1]['canonicalName']
		gbif_taxa += g_parents + [g_org]
		if len(gbif_taxa) > 250 and not _th.is_alive():
			_th = Thread(target=push_taxonomy, name='pushing taxonomy', args=['gbif', gbif_taxa])
			_th.start()
			gbif_taxa = []

	push_taxonomy('biota', gbif_taxa)
	return 'DONE'
