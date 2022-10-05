import json
from threading import Thread

from .. import create_request, read_request
from ....io.taxa import TaxaTaskTools, BIOTA_COLUMNS


def run():

	print("Initializing organism entries with Biota taxa")
	df = TaxaTaskTools.biota_get_df()
	df = df[df.kingdom == 'Plantae'][df.environment == 'Terrestre']		# only the Plantae kingdom (species: 3791 instead of 25161)

	print(' > Getting Biota species')
	for i, org in df[~df.species.str.contains('ssp.') | ~df.species.str.contains('subsp.')].iterrows():
		print(create_request('organisms/', split_name=1, **org, rank='species'))
	_th = Thread(target=read_request, name='Look for canonical names',
				args=['organisms/'], kwargs={'filter': {'rank': 'species'}})
	_th.start()

	print(' > Getting Biota subspecies')
	for i, org in df[df.species.str.contains('ssp.') | df.species.str.contains('subsp.')].iterrows():
		print(create_request('organisms/', split_name=1, **org, rank='subspecies'))
	_th = Thread(target=read_request, name='Look for canonical names',
				args=['organisms/'], kwargs={'filter': {'rank': 'subspecies'}})
	_th.start()

	ranks = read_request('ontologies/terms/', filter={'cv': 'taxonomic_rank'})
	for rank in [_.get('name') for _ in json.loads(ranks.text).get('content')
					if _.get('name') and _['name'] != 'species' and _['name'] in BIOTA_COLUMNS.values()]:
		print(' > Getting Biota ' + rank)
		for org in TaxaTaskTools.biota_get_by_rank(df, rank):
			print(create_request('organisms/', species=org, rank=rank))
		_th = Thread(target=read_request, name='Look for canonical names',
					args=['organisms/'], kwargs={'filter': {'rank': rank}})
		_th.start()

	return 'DONE'
