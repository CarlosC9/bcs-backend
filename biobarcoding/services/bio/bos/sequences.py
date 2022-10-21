import re

from Bio import SeqIO

from . import BosService
from ..meta.ontologies import get_type_id
from ..meta.organisms import Service as OrgService
from ... import get_orm_params, get_or_create, force_underscored, log_exception
from ...main import get_orm
from ....db_models import DBSession
from ....db_models.chado import Organism, StockFeature
from ....db_models.bioinformatics import Sequence, Specimen
from ....rest import filter_parse

org_service = OrgService()


##
# SEQUENCE SERVICE
##
class Service(BosService):

    def __init__(self):
        super(Service, self).__init__()
        self.orm = get_orm('sequences')
        self.obj_type = 'sequence'
        self.fos = Sequence
        self.formats = ['clustal', 'embl', 'fasta', 'genbank', 'gb', 'nexus', 'phylip', 'seqxml', 'abi', 'abi-trim', 'ace', 'cif-atom', 'cif-seqres', 'fasta-2line', 'fastq-sanger', 'fastq', 'fastq-solexa', 'fastq-illumina', 'gck', 'ig', 'imgt', 'pdb-seqres', 'pdb-atom', 'phd', 'pir', 'sff', 'sff-trim', 'snapgene', 'stockholm', 'swiss', 'tab', 'qual', 'uniprot-xml', 'xdna']
        # only write in 'clustal', 'embl', 'fasta', 'fasta-2line', 'fastq-sanger', 'fastq', 'fastq-solexa', 'fastq-illumina', 'genbank', 'gb', 'imgt', 'nexus', 'phd', 'phylip', 'pir', 'seqxml', 'sff', 'stockholm', 'tab', 'qual', 'xdna'

    def prepare_values(self, **values):

        if values.get('residues') and not values.get('seqlen'):
            values['seqlen'] = len(values.get('residues'))

        if not values.get('organism_id') and values.get('organism'):
            try:
                orgs = org_service.read(filter={'name': values.get('organism')})[0]
                values['organism_id'] = orgs[0].organism_id
            except Exception as e:
                values['organism_id'] = org_service.create(organism=values.get('organism'), split_name=1)[0].organism_id

        if not values.get('type_id') and values.get('type'):
            try:
                values['type_id'] = get_type_id(type=values.get('type'), subtype=values.get('subtype'))
            except:
                pass

        return super(Service, self).prepare_values(**values)

    def check_values(self, **values):

        if not values.get('uniquename'):
            raise Exception('Missing the uniquename (sequences ID)')

        if not values.get('organism_id'):
            values['organism_id'] = get_or_create(self.db, Organism, genus='unknown', species='organism').organism_id

        if not values.get('type_id'):
            try:
                values['type_id'] = get_type_id(type='sequence')
            except:
                raise Exception(f'Missing the type_id for {values.get("uniquename")}')

        return get_orm_params(self.orm, **values)

    def seq_stock(self, seq, **values):
        if not values.get('stock_id') and values.get('stock'):
            from ..meta.individuals import Service as StockService
            try:
                stock = StockService().get_query(purpose='annotate',
                                                 id=values.get('stock_id'),
                                                 uniquename=values.get('stock'))[0].one()
            except:
                stock = StockService().create(uniquename=values.get('stock'), organism_id=seq.organism_id)[0]
                self.db.flush()   # stock_id required
            values['stock_id'] = stock.stock_id
        elif not values.get('stock_id'):
            return None
        stock_bind = get_or_create(self.db, StockFeature,
                                   stock_id=values.get('stock_id'),
                                   feature_id=seq.feature_id,
                                   type_id=seq.type_id)
        return stock_bind

    def after_create(self, new_object, **values):
        values = super(Service, self).after_create(new_object, **values)

        stock = self.seq_stock(new_object, **values)
        # seq to bcs
        fos = get_or_create(DBSession, self.fos,
                            specimen_id=None if not stock else
                                DBSession.query(Specimen.id).filter(Specimen.native_id == stock.stock_id).one(),
                            native_id=new_object.feature_id,
                            name=new_object.uniquename)

        return values

    ##
    # READ
    ##

    def attach_data(self, *content) -> list:
        new = super(Service, self).attach_data(*content)

        _ids = [_['feature_id'] for _ in new]
        seqs = {}
        try:
            # from sqlalchemy.sql.expression import case
            # _ord = case({_id: index for index, _id in enumerate(_ids)},
            #             value=self.fos.native_id)   # TODO test order
            seqs = dict(DBSession.query(self.fos.native_id, self.fos.uuid)
                        .filter(self.fos.native_id.in_(_ids)).all())
        except Exception as e:
            print('Error: Additional data could not be attached.')
            log_exception(e)
        for _ in new:
            _['uuid'] = str(seqs.get(_['feature_id'], ''))

        return new

    ##
    # DELETE
    ##

    def delete_related(self, *content, **kwargs):
        ids = [seq.feature_id for seq in content]
        query = DBSession.query(self.fos).filter(self.fos.native_id.in_(ids))
        return len([DBSession.delete(row) for row in query.all()])

    ##
    # DEPRECATED IMPORT
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
            features = [{'type': f.split()[-1], 'qualifiers': {f.split()[-1]:f[:-len(f.split()[-1])].strip()}} for f in features]
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

    ##
    # IMPORT
    ##

    def read_infile(self, file, _format) -> any:
        return next(SeqIO.parse(file, _format))

    def import_file(self, infile, format=None, **kwargs):
        content_file, _format = self.check_infile(infile, format)
        from ....io.sequences import import_file
        return import_file(infile, _format, **kwargs)

    ##
    # EXPORT
    ##

    def seqs_header_parser(self, seqs, format: str):     # return dict(uniquename, header)
        # TODO: seqs header parser
        headers = {}
        if format:
            if "organism" in format.lower():    # ('organismID', 'organism', 'organism_canon', 'organism_canon_underscored')
                orgs = self.db.query(self.orm.uniquename, Organism.organism_id, Organism.name) \
                    .join(Organism).filter(self.orm.uniquename.in_([x.uniquename for x in seqs])).all()
                from ...species_names import parse_species_names
                underscored = "underscore" in format.lower()
                for seq_id, org_id, name in orgs:
                    if "id" in format.lower():
                        headers[seq_id] = str(org_id)
                    elif "canon" in format.lower():
                        headers[seq_id] = parse_species_names(DBSession, [name],
                                                              underscores=underscored)[0]
                    if not headers.get(seq_id):
                        headers[seq_id] = force_underscored(name) if underscored else name
        return headers

    def chado2biopy(self, seqs: list, header: str = None) -> list:
        headers = self.seqs_header_parser(seqs, header)
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        records = []
        for seq in seqs:
            # TODO: study molecule_type ?, and append taxonomy and features
            annotations = {
                'molecule_type': 'DNA',
                'organism': org_service.get_org_name(seq.organism_id)
            }
            records.append(SeqRecord(Seq(seq.residues or ''),
                                     headers.get(seq.uniquename, seq.uniquename),
                                     seq.uniquename,
                                     '' if headers else seq.name or '',
                                     annotations=annotations))
        return records

    def data2file(self, data: list, outfile, format: str, values={}, **kwargs) -> int:
        records = self.chado2biopy(data, values.pop('header', ''))
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
        return len(data)

    ##
    # GETTER AND OTHERS
    ##

    def aux_filter(self, _filter: dict) -> list:
        clauses = []

        if _filter.get('genus'):
            from ....db_models.chado import Organism
            _ids = self.db.query(Organism.organism_id)\
                .filter(filter_parse(Organism, [{'genus': _filter.get('genus')}]))
            clauses.append(self.orm.organism_id.in_(_ids))

        if _filter.get('stock_id'):
            from ....db_models.chado import StockFeature
            _ids = self.db.query(StockFeature.feature_id)\
                .filter(filter_parse(StockFeature, [{'stock_id': _filter.get('stock_id')}]))
            clauses.append(self.orm.feature_id.in_(_ids))

        if _filter.get('analysis_id'):
            from ....db_models.chado import AnalysisFeature
            _ids = self.db.query(AnalysisFeature.feature_id)\
                .filter(filter_parse(AnalysisFeature, [{'analysis_id': _filter.get('analysis_id')}]))
            clauses.append(self.orm.feature_id.in_(_ids))

        if _filter.get('phylotree_id'):
            from ....db_models.chado import Phylonode
            _ids = self.db.query(Phylonode.feature_id)\
                .filter(filter_parse(Phylonode, [{'phylotree_id': _filter.get('phylotree_id')}]))
            clauses.append(self.orm.feature_id.in_(_ids))

        if _filter.get("prop_cvterm_id"):
            from ....db_models.chado import Featureprop
            _ids = self.db.query(Featureprop.feature_id)\
                .filter(filter_parse(Featureprop, [{'type_id': _filter.get('prop_cvterm_id')}]))
            clauses.append(self.orm.feature_id.in_(_ids))

        if _filter.get("program"):
            from ....db_models.chado import Analysis, AnalysisFeature
            _ids = self.db.query(AnalysisFeature.feature_id).join(Analysis) \
                .filter(filter_parse(Analysis, [{'program': _filter.get('program')}]))
            clauses.append(self.orm.feature_id.in_(_ids))

        if _filter.get("programversion"):
            from ....db_models.chado import Analysis, AnalysisFeature
            _ids = self.db.query(AnalysisFeature.feature_id).join(Analysis) \
                .filter(filter_parse(Analysis, [{'programversion': _filter.get('programversion')}]))
            clauses.append(self.orm.feature_id.in_(_ids))

        if _filter.get("algorithm"):
            from ....db_models.chado import Analysis, AnalysisFeature
            _ids = self.db.query(AnalysisFeature.feature_id).join(Analysis) \
                .filter(filter_parse(Analysis, [{'algorithm': _filter.get('algorithm')}]))
            clauses.append(self.orm.feature_id.in_(_ids))

        from datetime import datetime
        if _filter.get("added-from"):
            _ids = self.db.query(self.orm.feature_id) \
                .filter(filter_parse(self.orm, {'timeaccessioned': _filter.get("added-from")}))
            clauses.append(self.orm.feature_id.in_(_ids))
        if _filter.get("added-to"):
            _ids = self.db.query(self.orm.feature_id) \
                .filter(filter_parse(self.orm, {'timeaccessioned': _filter.get("added-to")}))
            clauses.append(self.orm.feature_id.in_(_ids))
        if _filter.get("lastmodified-from"):
            _ids = self.db.query(self.orm.feature_id) \
                .filter(filter_parse(self.orm, {'timelastmodified': _filter.get("lastmodified-from")}))
            clauses.append(self.orm.feature_id.in_(_ids))
        if _filter.get("lastmodified-to"):
            _ids = self.db.query(self.orm.feature_id) \
                .filter(filter_parse(self.orm, {'timelastmodified': _filter.get("lastmodified-to")}))
            clauses.append(self.orm.feature_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(_filter)
