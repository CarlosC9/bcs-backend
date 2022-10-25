import json
import pandas as pd
from redis.client import Redis
from redis_lock import Lock

from biobarcoding import app_acronym
from .. import create_request, read_request, update_request, get_response_content, get_response_count, get_response_id
from ... import get_redis_host

NAMES = {
    'gbif': 'GBIF taxonomy tree',
    'ncbi': 'NCBI taxonomy tree',
    'biota': 'Biota taxonomy tree',
}

RANKS = ('superkingdom', 'kingdom', 'superphylum', 'phylum', 'subphylum', 'superclass', 'class', 'subclass',
         'superorder', 'order', 'suborder', 'superfamily', 'family', 'subfamily', 'genus', 'subgenus',
         'species', 'subspecies')

TAXA_LOCK = Lock(Redis(get_redis_host()), f"{app_acronym}-sysadmin-lock")


def get_taxonomic_ranks() -> any:
	ranks = read_request('ontologies/terms/', filter={'cv': ('taxonomy', 'taxonomic_rank')})
	return json.loads(ranks.text).get('content')


def get_ngd_genus() -> set:
	try:
		TAXA_LOCK.acquire()
		ngd_taxa = [t for t in get_response_content(
			read_request('organisms/', filter={'genus': {'op': 'ne', 'unary': ''}}, no_attach=True))]
		return set(pd.DataFrame(ngd_taxa)['genus'])
	finally:
		if TAXA_LOCK.locked():
			TAXA_LOCK.release()


def push_taxonomy(taxonomy: str = 'gbif', taxa: list = []):

	print('push_taxonomy:', len(taxa))
	g_df = pd.DataFrame(taxa)[['canonicalName', 'parent', 'rank']].drop_duplicates()
	g_df = g_df.where(pd.notnull(g_df), None)
	top_data = g_df[g_df.parent.isnull()]		# .to_dict('records')

	def l_recursive(i_data):
		for i, r in i_data.iterrows():
			values = {'species': r['canonicalName'], 'rank': r['rank'].lower()}
			params = {'taxonomy': NAMES.get(taxonomy, taxonomy), 'parent': r['parent']}
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
