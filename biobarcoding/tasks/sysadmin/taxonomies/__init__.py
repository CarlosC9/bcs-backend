import json

from .. import read_request


def get_taxonomic_ranks():
	ranks = read_request('ontologies/terms/', filter={'cv': ('taxonomy', 'taxonomic_rank')})
	return json.loads(ranks.text).get('content')
