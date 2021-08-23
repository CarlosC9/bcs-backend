from . import log_exception
from ..db_models import DBSessionChado as chado_session
from ..db_models.chado import Cvterm, Cv
from ..rest import Issue, IType, get_query


##
# CVTERMS: STATIC GETTERS
##

def get_seq_cvterm(type='sequence'):
    return {
        'sequence':                 {'cv':'sequence', 'name':'sequence_feature'},
        'partial sequence':         {'cv':'sequence', 'name':'region'},
        'gene':                     {'cv':'sequence', 'name':'gene'},
        'partial gene':             {'cv':'sequence', 'name':'gene_fragment'},
        'cds':                      {'cv':'sequence', 'name':'CDS'},
        'partial cds':              {'cv':'sequence', 'name':'CDS_fragment'},
        'consensus':                {'cv':'sequence', 'name':'consensus'},
        'processed':                {'cv':'sequence', 'name':'experimental_feature'},
        'aligned':                  {'cv':'sequence', 'name':'sequence_assembly'},
        'chloroplast':              {'cv':'sequence', 'name':'chloroplast_DNA'},
        'maturase k (matk) gene':   {'cv':'sequence', 'name':'gene', 'value':'matk'},
        'maturase k (rbcl) gene':   {'cv':'sequence', 'name':'gene', 'value':'rbcl'},
    }.get(type)


def get_aln_cvterm(type='alignment'):
    return {
        'alignment':    {'cv':'data', 'name':'Sequence alignment'},
        'hybrid':       {'cv':'data', 'name':'Hybrid sequence alignment'},
        'nucleic':      {'cv':'data', 'name':'Nucleic acid sequence alignment'},
        'pair':         {'cv':'data', 'name':'Pair sequence alignment'},
        'protein':      {'cv':'data', 'name':'Protein sequence alignment'},
    }.get(type)


def get_phy_cvterm(type='phylotree'):
    return {
        'phylotree':    {'cv':'data', 'name':'Phylogenetic tree'},
        'gene':         {'cv':'data', 'name':'Gene tree'},
        'species':      {'cv':'data', 'name':'Species tree'},
        'taxonomy':     {'cv':'data', 'name':'Taxonomy'},
    }.get(type)


def get_pda_cvterm(type='pda'):
    return {
        'pda':      {'cv':'data', 'name':'Biodiversity data'},      # Machine-readable biodiversity data.
        'alpha':    {'cv':'data', 'name':'Alpha diversity data'},   # The mean species diversity in sites or habitats at a local scale.
        'beta':     {'cv':'data', 'name':'Beta diversity data'},    # The ratio between regional and local species diversity.
        'gamma':    {'cv':'data', 'name':'Gamma diversity data'},   # The total species diversity in a landscape.
    }.get(type)


def get_blast_cvterm(type='blast'):
    return {
        'blast':    {'cv':'data', 'name':'Sequence search results'},
        'format':   {'cv':'format', 'name':'BLAST results'},
    }.get(type)


def get_tax_cvterm(type='taxonomy'):
    return {
        'taxonomy':     {'cv':'data', 'name':'Taxonomy'},
        'rank':         {'cv':'data', 'name':'Taxonomic classification'},
        'taxon':        {'cv':'data', 'name':'Taxon'},
        'kingdom':      {'cv':'data', 'name':'Kingdom name'},
        'family':       {'cv':'data', 'name':'Family name'},
        'genus':        {'cv':'data', 'name':'Genus name'},
        'species':      {'cv':'data', 'name':'Species name'},
        'ncbi':         {'cv':'data', 'name':'NCBI taxon'},
    }.get(type)


def get_ont_cvterm(type='ontology'):
    return {
        'ontology':     {'cv':'data', 'name':'Ontology'},
        'cv':           {'cv':'data', 'name':'Ontology'},
        'cvterm':       {'cv':'data', 'name':'Ontology concept'},
    }.get(type)


def get_used_cvterm(type, subtype=None):
    return {
        # 'default':      {'cv':'data', 'name':'Data'},
        'sequence':     get_seq_cvterm(subtype or type),
        # 'analysis':     {'cv':'', 'name':''},,
        'alignment':    get_aln_cvterm(subtype or type),
        'phylotree':    get_phy_cvterm(subtype or type),
        'pda':          get_pda_cvterm(subtype or type),
        'blast':        get_blast_cvterm(subtype or type),
        'taxonomy':     get_tax_cvterm(subtype or type),
        'ontology':     get_ont_cvterm(subtype or type),
    }.get(type)


##
# CV: CREATE
##

def __check_cv_values(**values):
    return dict([(k, v) for k, v in values.items() if k in Cv.__table__.columns])


def create(**kwargs):
    content = None
    try:
        values = __check_cv_values(**kwargs)
        chado_session.add(Cv(**values))
        issues, status = [Issue(IType.INFO, f'CREATE ontologies: The cv "{kwargs.get("name")}" was created successfully.')], 201
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, f'CREATE ontologies: The cv "{kwargs.get("name")}" could not be created.')], 409
    return issues, content, status


##
# CV: READ
##

def read(id=None, **kwargs):
    content = None
    try:
        content, count = get_cv_query(id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ ontologies: The ontologies were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ ontologies: The ontologies could not be read.')], 400
    return issues, content, count, status


##
# CV: UPDATE
##

def update(id, remote_url=None, input_file=None, **kwargs):
    content = None
    try:
        content = get_cv_query(id)[0].one()
        content.update(__check_cv_values(**kwargs))
        issues, status = [Issue(IType.INFO, f'UPDATE ontologies: The cv "{content.name}" was successfully updated.')], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, f'UPDATE ontologies: The cv "{content.name}" could not be updated.')], 409
    return issues, content, status


##
# CV: DELETE
##

def delete(id=None, **kwargs):
    content = None
    try:
        query, count = get_cv_query(id, **kwargs)
        content = query.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE ontologies: {count} cvs were successfully removed.')], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, f'DELETE ontologies: The cvs could not be removed.')], 404
    return issues, content, status


##
# CV: IMPORT
##

def import_file(input_file, format='obo', **kwargs):
    content = None
    try:
        from flask import current_app
        cfg = current_app.config
        # f"""go2fmt.pl -p obo_text -w xml {input_file} | \
        #     go-apply-xslt oboxml_to_chadoxml - > {input_file}.xml"""
        # cmd = f"""go2chadoxml {input_file} > /tmp/{input_file}.chado.xml;
        #     stag-storenode.pl -d 'dbi:Pg:dbname={cfg['database']};host={cfg['host']};port=cfg['port']'
        #     --user {cfg['user']} --password {cfg['password']} /tmp/{input_file}.chado.xml"""
        import pronto
        import os
        namespace = pronto.Ontology(input_file).metadata.default_namespace
        onto_name = f'-c {namespace}' if namespace else f'-c {os.path.basename(input_file)}'
        from biobarcoding.services import exec_cmds
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
        issues, status = [Issue(IType.INFO,
                                f'IMPORT ontologies: The {format} ontology was successfully imported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'IMPORT ontologies: The file {input_file} could not be imported.')], 409
    return issues, content, status


##
# CV: EXPORT
##

def export(id=None, format='obo', **kwargs):
    issues = [Issue(IType.WARNING, 'EXPORT ontologies: dummy completed')]
    return issues, None, 200


##
# CV: GETTER AND OTHERS
##

def get_cv_query(id=None, **kwargs):
    if id:
        query = chado_session.query(Cv).filter(Cv.cv_id == id)
        return query, query.count()
    return get_query(chado_session, Cv, aux_filter=__aux_own_filter, aux_order=__aux_own_order, **kwargs)


def __aux_own_filter(filter):
    clauses = []
    return clauses


def __aux_own_order(order):
    clauses = []
    return clauses


##
# CVTERMS: READ
##

def read_cvterms(cv_id=None, cvterm_id=None, **kwargs):
    content = None
    try:
        content, count = get_cvterm_query(cv_id=cv_id, id=cvterm_id, **kwargs)
        if cvterm_id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.ERROR, f'READ ontology terms: The cvterms were successfully read.')], 200
    except Exception as e:
        issues, status = [Issue(IType.ERROR, f'READ ontology terms: The cvterms could not be read.')], 400
    return issues, content, count, status


def get_cvterm_query(id=None, type=None, subtype=None, **kwargs):
    if id:
        query = chado_session.query(Cvterm).filter(Cvterm.cvterm_id == id)
        return query, query.count()
    if type:
        kwargs.update(**get_used_cvterm(type, subtype))
    if kwargs.get('cv'):
        kwargs['cv_id'] = chado_session.query(Cv.cv_id).filter(Cv.name==kwargs.get('cv'))
    return get_query(chado_session, Cvterm, aux_filter=__aux_cvterms_filter, aux_order=__aux_cvterms_order, **kwargs)


# TODO: featureprop, analysisprop, phylotreeprop, etc ?
def __aux_cvterms_filter(filter):
    clause = []
    if filter.get('feature_id'):
        from biobarcoding.db_models.chado import Feature, Featureprop, FeatureCvterm
        _ids = chado_session.query(Feature.type_id) \
            .filter(Feature.feature_id == filter.get('feature_id')).all()
        _ids += chado_session.query(Featureprop.type_id) \
            .filter(Featureprop.feature_id == filter.get('feature_id')).all()
        _ids += chado_session.query(FeatureCvterm.cvterm_id) \
            .filter(FeatureCvterm.feature_id == filter.get('feature_id')).all()
        clause.append(Cvterm.cvterm_id.in_(_ids))
    if filter.get('analysis_id'):
        from biobarcoding.db_models.chado import Analysisprop, AnalysisCvterm
        _ids = chado_session.query(Analysisprop.type_id) \
            .filter(Analysisprop.analysis_id == filter.get('analysis_id')).all()
        _ids += chado_session.query(AnalysisCvterm.cvterm_id) \
            .filter(AnalysisCvterm.analysis_id == filter.get('analysis_id')).all()
        clause.append(Cvterm.cvterm_id.in_(_ids))
    if filter.get('phylotree_id'):
        from biobarcoding.db_models.chado import Phylotree, Phylotreeprop
        _ids = chado_session.query(Phylotree.type_id) \
            .filter(Phylotree.phylotree_id == filter.get('phylotree_id')).all()
        _ids += chado_session.query(Phylotreeprop.type_id) \
            .filter(Phylotreeprop.phylotree_id == filter.get('phylotree_id')).all()
        clause.append(Cvterm.cvterm_id.in_(_ids))
    return clause


def __aux_cvterms_order(order):
    clauses = []
    return clauses
