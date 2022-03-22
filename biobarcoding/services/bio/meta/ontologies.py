from . import MetaService
from ... import log_exception
from ...main import get_orm
from ....db_models import DBSessionChado


##
# CVTERMS: STATIC GETTERS
##

def get_seq_cvterm(type='sequence'):
    return {
        'sequence':                 {'cv': 'sequence', 'name': 'sequence_feature'},
        'partial sequence':         {'cv': 'sequence', 'name': 'region'},
        'gene':                     {'cv': 'sequence', 'name': 'gene'},
        'partial gene':             {'cv': 'sequence', 'name': 'gene_fragment'},
        'cds':                      {'cv': 'sequence', 'name': 'CDS'},
        'partial cds':              {'cv': 'sequence', 'name': 'CDS_fragment'},
        'consensus':                {'cv': 'sequence', 'name': 'consensus'},
        'processed':                {'cv': 'sequence', 'name': 'experimental_feature'},
        'aligned':                  {'cv': 'sequence', 'name': 'sequence_assembly'},
        'chloroplast':              {'cv': 'sequence', 'name': 'chloroplast_DNA'},
        'matk':                     {'cv': 'sequence', 'name': 'gene', 'value': 'matk'},
        'maturase k (matk) gene':   {'cv': 'sequence', 'name': 'gene', 'value': 'matk'},
        'rbcl':                     {'cv': 'sequence', 'name': 'gene', 'value': 'rbcl'},
        'ribulose (rbcl) gene':     {'cv': 'sequence', 'name': 'gene', 'value': 'rbcl'},     # Ribulose bisphosphate carboxylase large
    }.get(type)


def get_stock_cvterm(type='stock'):
    return {
        'stock':        {'cv': 'sequence', 'name': 'wild_type'},
        'wild':         {'cv': 'sequence', 'name': 'wild_type'},
        'collection':   {'cv': 'sequence', 'name': 'variant_collection'},
        'genome':       {'cv': 'sequence', 'name': 'variant_genome'},
    }.get(type)


def get_stockcoll_cvterm(type='stockcoll'):
    return {
        'stockcoll':    {'cv': 'sequence', 'name': 'sequence_collection'},
        'sequence':     {'cv': 'sequence', 'name': 'sequence_collection'},
        'contig':       {'cv': 'sequence', 'name': 'contig_collection'},
    }.get(type)


def get_aln_cvterm(type='alignment'):
    return {
        'alignment':    {'cv': 'data', 'name': 'Sequence alignment'},
        'hybrid':       {'cv': 'data', 'name': 'Hybrid sequence alignment'},
        'nucleic':      {'cv': 'data', 'name': 'Nucleic acid sequence alignment'},
        'pair':         {'cv': 'data', 'name': 'Pair sequence alignment'},
        'protein':      {'cv': 'data', 'name': 'Protein sequence alignment'},
    }.get(type)


def get_phy_cvterm(type='phylotree'):
    return {
        'phylotree':    {'cv': 'data', 'name': 'Phylogenetic tree'},
        'gene':         {'cv': 'data', 'name': 'Gene tree'},
        'species':      {'cv': 'data', 'name': 'Species tree'},
        'taxonomy':     {'cv': 'data', 'name': 'Taxonomy'},
    }.get(type)


def get_pda_cvterm(type='pda'):
    return {
        'pda':      {'cv': 'data', 'name': 'Biodiversity data'},      # Machine-readable biodiversity data.
        'alpha':    {'cv': 'data', 'name': 'Alpha diversity data'},   # The mean species diversity in sites or habitats at a local scale.
        'beta':     {'cv': 'data', 'name': 'Beta diversity data'},    # The ratio between regional and local species diversity.
        'gamma':    {'cv': 'data', 'name': 'Gamma diversity data'},   # The total species diversity in a landscape.
    }.get(type)


def get_blast_cvterm(type='blast'):
    return {
        'blast':    {'cv': 'data', 'name': 'Sequence search results'},
        'format':   {'cv': 'format', 'name': 'BLAST results'},
    }.get(type)


def get_tax_cvterm(type='taxonomy'):
    return {
        'taxonomy':     {'cv': 'data', 'name': 'Taxonomy'},
        'rank':         {'cv': 'data', 'name': 'Taxonomic classification'},
        'taxon':        {'cv': 'data', 'name': 'Taxon'},
        'kingdom':      {'cv': 'data', 'name': 'Kingdom name'},
        'family':       {'cv': 'data', 'name': 'Family name'},
        'genus':        {'cv': 'data', 'name': 'Genus name'},
        'species':      {'cv': 'data', 'name': 'Species name'},
        'ncbi':         {'cv': 'data', 'name': 'NCBI taxon'},
    }.get(type)


def get_ont_cvterm(type='ontology'):
    return {
        'ontology':     {'cv': 'data', 'name': 'Ontology'},
        'cv':           {'cv': 'data', 'name': 'Ontology'},
        'cvterm':       {'cv': 'data', 'name': 'Ontology concept'},
    }.get(type)


def get_used_cvterm(type, subtype=None):
    return {
        # 'default':      {'cv':'data', 'name':'Data'},
        'sequence':     get_seq_cvterm(subtype or type),
        'stock':        get_stock_cvterm(subtype or type),
        'stockcoll':    get_stockcoll_cvterm(subtype or type),
        # 'analysis':     {'cv':'', 'name':''},,
        'alignment':    get_aln_cvterm(subtype or type),
        'phylotree':    get_phy_cvterm(subtype or type),
        'pda':          get_pda_cvterm(subtype or type),
        'blast':        get_blast_cvterm(subtype or type),
        'taxonomy':     get_tax_cvterm(subtype or type),
        'ontology':     get_ont_cvterm(subtype or type),
    }.get(type)


def get_type_id(type=None, subtype=None):
    return CvtermService().get_query(**get_used_cvterm(type=type, subtype=subtype))[0].one().cvterm_id


##
# ONTOLOGY SERVICE
##
class Service(MetaService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('ontologies')

    ##
    # IMPORT
    ##

    def import_file(self, input_file, format='obo', **kwargs):
        content, count = None, 0
        try:
            from flask import current_app
            cfg = current_app.config
            # f"""go2fmt.pl -p obo_text -w xml {input_file} | \
            #     go-apply-xslt oboxml_to_chadoxml - > {input_file}.xml"""
            # cmd = f"""go2chadoxml {input_file} > /tmp/{input_file}.chado.xml;
            #     stag-storenode.pl -d 'dbi:Pg:dbname={cfg['database']};host={cfg['host']};port=cfg['port']'
            #     --user {cfg['user']} --password {cfg['password']} /tmp/{input_file}.chado.xml"""
            import pronto
            import os.path
            namespace = pronto.Ontology(input_file).metadata.default_namespace
            onto_name = f'-c {namespace}' if namespace else f'-c {os.path.basename(input_file)}'
            from ... import exec_cmds
            out, err = exec_cmds([
                f'''perl ./biobarcoding/services/perl_scripts/gmod_load_cvterms.pl\
                    -H {cfg["CHADO_HOST"]}\
                    -D {cfg["CHADO_DATABASE"]}\
                    -r {cfg["CHADO_USER"]}\
                    -p {cfg["CHADO_PASSWORD"]}\
                    -d Pg -s null -u\
                    {input_file}''',
                f'''perl ./biobarcoding/services/perl_scripts/gmod_make_cvtermpath.pl\
                    -H {cfg["CHADO_HOST"]}\
                    -D {cfg["CHADO_DATABASE"]}\
                    -u {cfg["CHADO_USER"]}\
                    -p {cfg["CHADO_PASSWORD"]}\
                    -d Pg {onto_name}'''])
            count += 1
        except Exception as e:
            log_exception(e)
            raise Exception(f'IMPORT ontologies: The file {os.path.basename(input_file)} could not be imported.')
        return content, count


##
# CVTERM SERVICE
##
class CvtermService(MetaService):

    def __init__(self):
        super(CvtermService, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('cvterms')

    def prepare_values(self, type=None, subtype=None, **values):

        if not values.get('cv_id') and values.get('cv'):
            from ....db_models.chado import Cv
            values['cv_id'] = self.db.query(Cv.cv_id).filter(Cv.name == values.get('cv'))

        if type:
            values.update(**get_used_cvterm(type, subtype))

        return super(CvtermService, self).prepare_values(**values)

    def aux_filter(self, filter):
        clauses = []

        if filter.get('feature_id'):
            from ....db_models.chado import Feature, Featureprop, FeatureCvterm
            _ids = self.db.query(Feature.type_id) \
                .filter(Feature.feature_id == filter.get('feature_id')).all()
            _ids += self.db.query(Featureprop.type_id) \
                .filter(Featureprop.feature_id == filter.get('feature_id')).all()
            _ids += self.db.query(FeatureCvterm.cvterm_id) \
                .filter(FeatureCvterm.feature_id == filter.get('feature_id')).all()
            clauses.append(self.orm.cvterm_id.in_(_ids))
        
        if filter.get('analysis_id'):
            from ....db_models.chado import Analysisprop, AnalysisCvterm
            _ids = self.db.query(Analysisprop.type_id) \
                .filter(Analysisprop.analysis_id == filter.get('analysis_id')).all()
            _ids += self.db.query(AnalysisCvterm.cvterm_id) \
                .filter(AnalysisCvterm.analysis_id == filter.get('analysis_id')).all()
            clauses.append(self.orm.cvterm_id.in_(_ids))
        
        if filter.get('phylotree_id'):
            from ....db_models.chado import Phylotree, Phylotreeprop
            _ids = self.db.query(Phylotree.type_id) \
                .filter(Phylotree.phylotree_id == filter.get('phylotree_id')).all()
            _ids += self.db.query(Phylotreeprop.type_id) \
                .filter(Phylotreeprop.phylotree_id == filter.get('phylotree_id')).all()
            clauses.append(self.orm.cvterm_id.in_(_ids))

        return clauses + super(CvtermService, self).aux_filter(filter)
