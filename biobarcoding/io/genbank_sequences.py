from Bio import Entrez, SeqIO, SeqRecord, SeqFeature

from biobarcoding import app_acronym

Entrez.email = f"admin@{app_acronym}.eu"
TAXON_SLOT_SIZE = "10"


class GenbankSeqsTools:

	@staticmethod
	def get_gbk_seqs(org: str, genes: tuple = ('matk', 'rbcl', 'its'), outfile="gbk_seqs.gb"):
		"""
		Get sequences of nucleotides from Genbank that match with a taxon and the given genes
		:param org: taxon name
		:param genes: list of gene names
		:param outfile: pathfile to store the data
		:return: outfile
		"""
		query = org + "[ORGN]" if org else 'plants[filter]'
		if genes:
			query += " AND (%s)" % (" OR ".join([f"{gene}[GENE]" for gene in genes]))
		handle = Entrez.esearch(db="nucleotide", retmax=TAXON_SLOT_SIZE, term=query)
		id_list = Entrez.read(handle)['IdList']
		if id_list:
			with Entrez.efetch(db='nucleotide', id=id_list, rettype='gb') as net_hdl:
				with open(outfile, "a") as out_handle:
					out_handle.write(net_hdl.read())
			print(f"> {org} saved.")
		else:
			print(f"! {org} not found.")
		return outfile

	@staticmethod
	def split_seqs(src_file: str, dst_file: str, genes: list = ('matk', 'rbcl', 'its'), _format: str = 'gb'):
		"""
		Extract genes from Genbank sequence file to a new file
		:param src_file: pathfile with the sequences
		:param dst_file: pathfile to store the data
		:param genes: genes to be extracted
		:param _format: input file format
		:return: dst_file
		"""

		def is_gene(feat: SeqFeature):
			"""
			Check if a SeqFeature is a gene, and it is wanted
			:param feat: a SeqFeature
			:return: bool
			"""
			if not feat.qualifiers.get('gene'):
				return False        # TODO: what is it ?
			return bool(len([_ for _ in feat.qualifiers.get('gene') if _.lower() in genes]))

		def extract_gene(seq: SeqRecord, gene: SeqFeature):
			"""
			Do something
			:param :
			:return:
			"""
			_ = gene.extract(seq)
			_.annotations = seq.annotations.copy()
			_.dbxrefs = seq.dbxrefs[:]
			_.id = '%s.%s' % (seq.id, gene.qualifiers.get('gene'))
			_.name, _.description = seq.name, seq.description
			return _

		result, missing = {}, 0
		with open(src_file) as input_handle:
			for record in SeqIO.parse(input_handle, _format):
				_len = len(result)
				for f in record.features:
					if f.type == "gene" and is_gene(f):
						new = extract_gene(record, f)
						if new.id not in result:
							result[new.id] = new
						# if "is complement":        # TODO: is it needed ?
						#     new = new.reverse_complement()
				if _len == len(result):
					missing += 1
		print('EXTRACTED: %s\tMISSING: %s' % (len(result), missing))

		with open(dst_file, "w") as fs:
			SeqIO.write(result.values(), fs, _format)
		print('Sequences in file %s splited.' % (src_file))

		return dst_file


def run():

	print("Initializing sequence entries from Genbank by Biota taxa")

	srcfile = '/tmp/gbk_seqs.gb'
	open(srcfile, 'w').close()

	print(' > Getting Biota species')
	from biobarcoding.tasks.sysadmin.initialize.taxa import TaxaTaskTools
	df = TaxaTaskTools.biota_get_df()
	df = df[df.kingdom == 'Plantae'][df.environment == 'Terrestre']

	print(' > Getting Genbank sequences')
	for org in df['species'].values:
		GenbankSeqsTools.get_gbk_seqs(" ".join(org.strip().split()[:2]), outfile=srcfile)
		GenbankSeqsTools.get_gbk_seqs(org, outfile=srcfile)

	print(' > Getting extracting genes')
	dstfile = '/tmp/gbk_sub_seqs.gb'

	GenbankSeqsTools.split_seqs(srcfile, dstfile)

	# TODO: import file

	return 'DONE'
