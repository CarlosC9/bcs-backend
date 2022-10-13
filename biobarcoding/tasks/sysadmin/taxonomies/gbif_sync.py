from threading import Thread
from pygbif import species

import pandas as pd

# from biobarcoding.tasks.sysadmin.taxonomies import get_taxonomic_ranks
from biobarcoding.tasks.sysadmin import create_request, read_request, update_request, \
	get_response_count, get_response_content, get_response_id


def push_taxonomy(*taxa):

	print('push_taxonomy:', len(taxa))
	g_df = pd.DataFrame(taxa)[['canonicalName', 'parent', 'rank']].drop_duplicates()
	g_df = g_df.where(pd.notnull(g_df), None)
	top_data = g_df[g_df.parent.isnull()]		# .to_dict('records')

	def l_recursive(i_data):
		for i, r in i_data.iterrows():
			values = {'species': r['canonicalName'], 'rank': r['rank'].lower()}
			params = {'taxonomy': 'GBIF Taxonomy tree', 'parent': r['parent']}
			try:
				response = read_request('organisms/', filter=values, no_attach=True)
				if get_response_count(response) != 1:
					if response.status_code < 300:
						print('organism not found, creating', values)
						response = create_request('organisms/', values={**values, **params})
					else:
						raise Exception()		# TODO PUT /api/auth ?
				response = update_request('organisms/%s' % get_response_id(response, 'organism_id'), values=params)
			except Exception as e:
				print('missing auth token ?', values)
			l_recursive(g_df[g_df.parent == r['canonicalName']])

	l_recursive(top_data)


def run():

	print("Synchronizing taxonomic hierarchy with GBIF")
	print(">> Get data from ngd")
	# ranks = get_taxonomic_ranks()
	ngd_taxa = [t for t in get_response_content(
		read_request('organisms/', filter={'genus': {'op': 'ne', 'unary': ''}}, no_attach=True))]
	ngd_genus = set(pd.DataFrame(ngd_taxa)['genus'])

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
			_th = Thread(target=push_taxonomy, name='pushing taxonomy', args=gbif_taxa)
			_th.start()
			gbif_taxa = []

	return 'DONE'
