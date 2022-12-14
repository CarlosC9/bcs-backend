import re
import os.path
import pandas as pd
from Bio import Entrez
import requests
from pygbif import species

from biobarcoding import app_acronym
from biobarcoding.services import get_encoding
from biobarcoding.services.species_names import get_taxon_from_gbif, get_lineage_from_gbif, get_children_from_gbif

EMAIL = f"admin@{app_acronym}.eu"
BIOTA_URL = "https://www.biodiversidadcanarias.es/biota/especies/export"
# BIOTA_URL_OPTIONS = "?pagina=1&tipoBusqueda=NOMBRE&searchSpeciesTabs=fastSearchTab&orderBy=nombreCientifico&orderForm=true"
# "?pagina=1&codigo=&tipoBusqueda=NOMBRE&nombreCientifico=&nombreComun=&subnomine=&tipo=2&reinoPk=3&divisionPk=9&subdivisionPk=6&clasePk=38&ordenPk=262&familiaPk=896&generoPk=4517&searchSpeciesTabs=scientificTab&orderBy=nombreCientifico&orderForm=true"
NCBI_URL = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lin=s&name='
GBIF_URL = 'https://api.gbif.org/v1/species/'
GBIF_MATCH_URL = GBIF_URL + 'match?name='

BIOTA_COLUMNS = {
		"Código": "code",
		"Especie": "species",
		"Nombre común/vulgar": "common_name",
		"Medio": "environment",
		"Clases": "class",
		"Familia": "family",
		"Origen": "origin",
		"Endemicidad": "endemicity",
		"Imagen": "image",
		"Reino": "kingdom",
		"Clasificación Fungi": "fungi_classification",
		"División": "division",
		"Subdivisión": "subdivision",
		"Filo": "phylum",
		"Orden": "order",
		"Presencia": "presence",
		"Validez": "validity",
		"Status": "status",
		"Género": "genus",
		"Código EU-Nomen": "EU-Nomen_code",
		"Código EUNIS": "EUNIS_code",
		"Código Natura2000": "Natura2000_code",
	}


class TaxaTaskTools:

	@staticmethod
	def biota_get_csv(file="/tmp/biota_species.csv", overwrite=True) -> str:
		"""
		Download and store the Biota taxa file
		:param file: pathfile to store the data
		:param overwrite: flag to overwrite the file
		:return: file
		"""
		print(">> WORKING ON: biota_get_csv")
		if overwrite or not os.path.isfile(file):
			r = requests.get(BIOTA_URL, verify=False)
			# r.encoding = r.apparent_encoding
			with open(file, 'w') as f:
				f.write(r.content.decode(r.encoding))
		return file

	@staticmethod
	def biota_parse_columns(df):
		"""
		Translate a Biota Dataframe
		:param df: a Biota Dataframe
		:return: translation of a Biota Dataframe
		"""
		return [BIOTA_COLUMNS.get(c, c) for c in df.keys()]

	@staticmethod
	def biota_get_df():
		"""
		Get BIOTA entries as Dataframe
		:return: BIOTA entries as Dataframe
		"""
		# # TODO: figure out the encoding and the delimiter
		# # r.apparent_encoding, r.encoding, chardet.detect(r.raw)
		# encodings = ['Cp1252', 'Windows-1252', 'ISO-8859-9', 'ascii']
		# for encoding in encodings:
		# 	try:
		# 		frame = pd.read_csv(BIOTA_URL, encoding=encoding, sep=';')
		# 		frame.rename(columns=BIOTA_COLUMNS, inplace=True)
		# 		return frame
		# 	except Exception as e:
		# 		pass
		try:
			_file = TaxaTaskTools.biota_get_csv(overwrite=False)
			frame = pd.read_csv(_file, encoding=get_encoding(_file), sep=';')
			frame.rename(columns=BIOTA_COLUMNS, inplace=True)
			return frame.where(frame.notnull(), None)
		except Exception as e:
			return []

	@staticmethod
	def biota_get_species(df=None):
		"""
		Get the list of taxa
		:param df: a Biota Dataframe
		:return: a list of taxa
		"""
		if not df:
			df = TaxaTaskTools.biota_get_df()
		return df['species'].values

	@staticmethod
	def biota_get_by_rank(df, rank):
		"""
		Get the list of taxa
		:param df: a Biota Dataframe
		:param rank: the desired rank
		:return: a list of taxa
		"""
		if df is None:
			df = TaxaTaskTools.biota_get_df()
		for taxa in df[rank].unique():
			if taxa:
				yield taxa

	@staticmethod
	def gbif_get_df(*orgs):
		"""
		Get GBIF entries as Dataframe
		:param orgs: list of taxa
		:return: GBIF entries as Dataframe
		"""
		if not orgs:
			orgs = TaxaTaskTools.biota_get_species()
		return pd.DataFrame([get_taxon_from_gbif(org) for org in orgs])

	@staticmethod
	def gbif_get_parents_df(*orgs):
		"""
		Get GBIF entries as Dataframe
		:param orgs: list of taxa
		:return: GBIF entries as Dataframe
		"""
		if not orgs:
			df = TaxaTaskTools.biota_get_df()
		else:
			df = TaxaTaskTools.gbif_get_df(*orgs)
		lineages = []
		for i, r in df.iterrows():
			parents = get_lineage_from_gbif(r['usageKey'].astype(int))
			lineages += parents + [{**r.to_dict(), 'parent': parents[-1]['canonicalName']}]
		return pd.DataFrame(lineages).drop_duplicates()

	@staticmethod
	def ncbi_taxa2txids(*orgs) -> list:
		"""
		Get NCBI txids
		:param orgs: list of taxa
		:return: GBIF entries as Dataframe
		"""
		if not orgs:
			orgs = TaxaTaskTools.biota_get_species()
		print(">> WORKING ON: ncbi_taxa2txids")
		Entrez.email = EMAIL
		term = " OR ".join([f"{'+'.join(org.split(' ')[:2])}[Orgn]" for org in orgs])
		with Entrez.esearch(db="taxonomy", term=term, usehistory="y", retmax=20) as handle:
			ez_search = Entrez.read(handle)
		txids = []
		batch_size = 1000
		for start in range(0, int(ez_search['Count']), batch_size):
			with Entrez.efetch(db="taxonomy", retstart=start, retmax=batch_size,
							   webenv=ez_search["WebEnv"], query_key=ez_search["QueryKey"]) \
					as handle:
				ez_fetch = Entrez.read(handle)
			txids += ez_fetch
			print(len(txids))
		return txids

	@staticmethod
	def ncbi_taxa2txids_by_curl(*orgs) -> list:
		"""
		Get NCBI txids
		:param orgs: list of taxa
		:return: GBIF entries as Dataframe
		"""
		if not orgs:
			orgs = TaxaTaskTools.biota_get_species()
		print(">> WORKING ON: ncbi_taxa2txids_by_curl")
		txids = set()
		for org in orgs:
			url = f"{NCBI_URL}{'+'.join(org.split(' ')[:2])}"
			r = requests.get(url)
			txid = set(re.findall(r'\W(txid\d+)\D', str(r.content)))
			txids.update(txid)
			print(txid, len(txids), org)
		return list(txids)
