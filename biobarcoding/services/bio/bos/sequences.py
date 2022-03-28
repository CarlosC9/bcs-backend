import os.path
import re

from . import BosService
from ..meta.ontologies import get_type_id
from ... import get_orm_params, log_exception, get_bioformat, get_or_create
from ...main import get_orm
from ....db_models import DBSession
from ....db_models import DBSessionChado
from ....db_models.chado import Feature, Organism, StockFeature
from ....db_models.bioinformatics import Sequence
from ....rest import filter_parse


##
# SEQUENCE SERVICE
##
class Service(BosService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('sequences')
        self.bos = 'sequence'

    def prepare_values(self, **values):
        """
        Process the values 'organism', 'stock', and 'type' into valid IDs
        """

        if values.get('residues') and not values.get('seqlen'):
            values['seqlen'] = len(values.get('residues'))

        if not values.get('organism_id') and values.get('organism'):
            org = values.get('organism').strip()
            if org:
                # TODO: check canónical name (getting genus and species?)
                genus = org.split()[0]
                species = org[len(genus):].strip()
                values['organism_id'] = get_or_create(self.db, Organism, genus=genus, species=species).organism_id
        if not values.get('organism_id'):
            values['organism_id'] = get_or_create(self.db, Organism, genus='unknown', species='Unclassified organism').organism_id

        if not values.get('type_id') and values.get('type'):
            try:
                values['type_id'] = get_type_id(type=values.get('type'), subtype=values.get('subtype'))
            except:
                pass

        return super(Service, self).prepare_values(**values)

    def check_values(self, **values):
        """
        Fill in the empty not null fields whenever possible.
        """
        if not values.get('uniquename'):
            raise Exception('Missing the uniquename (sequences ID)')
        if not values.get('type_id'):
            try:
                values['type_id'] = get_type_id(type='sequence')
            except:
                raise Exception(f'Missing the type_id for {values.get("uniquename")}')
        return get_orm_params(self.orm, **values)

    def seq_stock(self, seq, **values):
        if not values.get('stock_id') and values.get('stock'):
            from ..meta.individuals import __get_query as get_stock_query, create as create_stock
            try:
                stock = get_stock_query(values.get('stock_id'), uniquename=values.get('stock'))[0].one()
            except:
                stock = create_stock(uniquename=values.get('stock'), organism_id=seq.organism_id)[1]
                self.db.flush()   # stock_id required
            values['stock_id'] = stock.stock_id
        elif not values.get('stock_id'):
            return None
        stock_bind = get_or_create(self.db, StockFeature,
                                   stock_id=values.get('stock_id'),
                                   feature_id=seq.feature_id,
                                   type_id=seq.type_id)
        return stock

    def after_create(self, new_object, **values):
        super(Service, self).after_create(new_object, **values)

        stock = self.seq_stock(new_object, **values)
        # seq to bcs
        fos_seq = get_or_create(DBSession, Sequence,   # specimen_id=fos_specimen.id
                                native_id=new_object.feature_id,
                                native_table='feature',
                                name=new_object.uniquename)
        return values

    ##
    # READ
    ##

    def attach_data(self, content):
        new = super(Service, self).attach_data(content)
        new['uuid'] = str(DBSession.query(Sequence)
                          .filter(Sequence.native_table == 'feature',
                                  Sequence.native_id == new.get('feature_id')).one().uuid)
        return new

    ##
    # DELETE
    ##

    def after_delete(self, *content, **kwargs):
        ids = [seq.feature_id for seq in content]
        query = DBSession.query(Sequence).filter(Sequence.native_id.in_(ids))
        return query.count()
        # return query.delete(synchronize_session='fetch')

    ##
    # IMPORT
    ##

    def simpleSeq2chado(self, seq, **params):
        return self.create(uniquename=seq.id, name=seq.description, residues=str(seq.seq), seqlen=len(seq.seq), **params)

    def fasSeq2chado(self, seq, **params):
        # params = { organism:<organism>, features=<features>, origin:<origin> }
        features = params.get('features')
        origin = params.get('origin')
        if features:
            if isinstance(features, str):
                features = features.split(',')
            features = [ {'type': f.split()[-1], 'qualifiers': {f.split()[-1]:f[:-len(f.split()[-1])].strip()}} for f in features ]
            params['features'] = features if isinstance(features, (tuple, list, set)) else params['features']
            params['molecule_type'] = features[-1]['type']
        elif origin:
            params['molecule_type'] = params.get('origin')
        return self.create(uniquename=seq.id, name=seq.description, stock=seq.name, residues=str(seq.seq), seqlen=len(seq.seq), **params)

    def gbSeq2chado(self, seq, **params):
        # seq.id, seq.name, seq.seq, seq.description, seq.dbxrefs, seq.letter_annotations
        # seq.annotations # (<class 'dict'>):
        # seq.features # (<class 'list'>):
        seq.name = f'{seq.id} {seq.description}' if seq.id == seq.name else seq.name
        return self.create(uniquename=seq.id, name=seq.name, residues=str(seq.seq), seqlen=len(seq.seq), description=seq.description,
                           dbxrefs=seq.dbxrefs, **seq.annotations, features=seq.features, **params)

    # ?META: accession, genus, species, molec_region, molec_origin, isle, georegion, individual, collection
    def bio2chado(self, seq, format, **params):
        try:
            if format == 'fasta':    # fasta
                pattern = '^(?P<id>\w+?) *\[organism=(?P<organism>.+?)\] *(?P<features>.+?); *(?P<origin>.+?)$'
                # p.e.: >Seq1_18037 [organism=Ruta montana] maturase K (matk) gene, partial sequence, partial cds; chloroplast
                meta = re.match(pattern, seq.description)
                if meta:
                    return self.fasSeq2chado(seq, **meta.groupdict(), **params)
            elif format == 'genbank':   # genbank
                return self.gbSeq2chado(seq, **params)
            return self.simpleSeq2chado(seq, **params)
        except:
            self.db.rollback()
            # return [Issue(IType.ERROR, f'IMPORT sequences: {seq.id} could not be imported.')]
            return None, 0

    def import_file(self, infile, format=None, **kwargs):
        content, count = [], 0
        format = get_bioformat(infile, format)
        try:
            from ... import seqs_parser
            for s in seqs_parser(infile, format):
                c, cc = self.bio2chado(s, format, **kwargs)
                content.append(c)
                count += cc
                # issues += [Issue(i.itype, i.message, os.path.basename(infile)) for i in iss]
        except Exception as e:
            log_exception(e)
            raise Exception(f'IMPORT sequences: file {os.path.basename(infile)} could not be imported.')
        return content, count

    ##
    # EXPORT
    ##

    def seqs_header_parser(self, seqs, format: str):     # return dict(uniquename, header)
        # TODO: seqs header parser
        headers = {}
        if "organism" in format.lower():    # ('organismID', 'organism', 'organism_canon', 'organism_canon_underscored')
            orgs = self.db.query(Feature.uniquename, Organism.organism_id, Organism.name) \
                .join(Organism).filter(Feature.uniquename.in_([x.uniquename for x in seqs])).all()
            from ...species_names import get_canonical_species_names
            for seq_id, org_id, name in orgs:
                if "id" in format.lower():
                    headers[seq_id] = str(org_id)
                elif "canon" in format.lower():
                    headers[seq_id] = get_canonical_species_names(DBSession, [name],
                                                    underscores="underscore" in format.lower())
                else:
                    headers[seq_id] = name
        return headers

    def chado2biopy(self, seqs: list, headers: dict) -> list:
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        records = []
        for seq in seqs:
            # TODO: study molecule_type ?, and append taxonomy and features
            annotations = {
                'molecule_type': 'DNA',
                'organism': self.db.query(Organism.name)
                    .filter(Organism.organism_id == seq.organism_id).one()[0]
            }
            records.append(SeqRecord(Seq(seq.residues),
                                     headers.get(seq.uniquename, seq.uniquename),
                                     seq.uniquename,
                                     '' if headers else seq.name,
                                     annotations=annotations))
        return records

    def data2file(self, seqs: list, outfile, format: str, values={}, **kwargs) -> int:
        headers = self.seqs_header_parser(seqs, values.pop('header', ''))
        from Bio import SeqIO
        records = self.chado2biopy(seqs, headers)
        res = SeqIO.write(records, outfile, format=format)
        if format == "nexus":
            # format datatype=dna missing=? gap=- matchchar=.;
            from Bio.Nexus.Nexus import Nexus
            aln = Nexus(outfile)
            aln.datatype = 'dna'
            aln.missing = '?'
            aln.gap = '-'
            aln.matchchar = '.'
            res = aln.write_nexus_data(outfile)
        return res

    ##
    # GETTER AND OTHERS
    ##

    def aux_filter(self, filter):
        clauses = []

        if 'analysis_id' in filter:
            from ....db_models.chado import AnalysisFeature
            _ids = self.db.query(AnalysisFeature.feature_id)\
                .filter(filter_parse(AnalysisFeature, [{'analysis_id': filter.get('analysis_id')}]))
            clauses.append(Feature.feature_id.in_(_ids))

        if 'phylotree_id' in filter:
            from ....db_models.chado import Phylonode
            _ids = self.db.query(Phylonode.feature_id)\
                .filter(filter_parse(Phylonode, [{'phylotree_id': filter.get('phylotree_id')}]))
            clauses.append(Feature.feature_id.in_(_ids))

        if "prop_cvterm_id" in filter:
            from ....db_models.chado import Featureprop
            _ids = self.db.query(Featureprop.feature_id)\
                .filter(filter_parse(Featureprop, [{'type_id': filter.get('prop_cvterm_id')}]))
            clauses.append(Feature.feature_id.in_(_ids))

        if "program" in filter:
            from ....db_models.chado import Analysis
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'program': filter.get('program')}]))
            from ....db_models.chado import AnalysisFeature
            _ids = self.db.query(AnalysisFeature.feature_id) \
                .filter(AnalysisFeature.analysis_id.in_(_ids))
            clauses.append(Feature.feature_id.in_(_ids))

        if "programversion" in filter:
            from ....db_models.chado import Analysis
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'programversion': filter.get('programversion')}]))
            from ....db_models.chado import AnalysisFeature
            _ids = self.db.query(AnalysisFeature.feature_id) \
                .filter(AnalysisFeature.analysis_id.in_(_ids))
            clauses.append(Feature.feature_id.in_(_ids))

        if "algorithm" in filter:
            from ....db_models.chado import Analysis
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'algorithm': filter.get('algorithm')}]))
            from ....db_models.chado import AnalysisFeature
            _ids = self.db.query(AnalysisFeature.feature_id) \
                .filter(AnalysisFeature.analysis_id.in_(_ids))
            clauses.append(Feature.feature_id.in_(_ids))

        from datetime import datetime
        if "added-from" in filter:
            filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
            _ids = self.db.query(Feature.feature_id) \
                .filter(filter_parse(Feature, {'timeaccessioned':filter.get("added-from")}))
            clauses.append(Feature.feature_id.in_(_ids))
        if "added-to" in filter:
            filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
            _ids = self.db.query(Feature.feature_id) \
                .filter(filter_parse(Feature, {'timeaccessioned':filter.get("added-to")}))
            clauses.append(Feature.feature_id.in_(_ids))
        if "lastmodified-from" in filter:
            filter["lastmodified-from"]['unary'] = datetime.strptime(filter.get("lastmodified-from")['unary'], '%Y-%m-%d')
            _ids = self.db.query(Feature.feature_id) \
                .filter(filter_parse(Feature, {'timelastmodified':filter.get("lastmodified-from")}))
            clauses.append(Feature.feature_id.in_(_ids))
        if "lastmodified-to" in filter:
            filter["lastmodified-to"]['unary'] = datetime.strptime(filter.get("lastmodified-to")['unary'], '%Y-%m-%d')
            _ids = self.db.query(Feature.feature_id) \
                .filter(filter_parse(Feature, {'timelastmodified':filter.get("lastmodified-to")}))
            clauses.append(Feature.feature_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(filter)
