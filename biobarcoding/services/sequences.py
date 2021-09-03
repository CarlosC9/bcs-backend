import os.path
import re

from .ontologies import get_cvterm_query
from ..db_models import DBSession as db_session
from ..db_models import DBSessionChado as chado_session
from ..db_models.chado import Feature, Organism, StockFeature
from ..db_models.bioinformatics import Sequence

from ..rest import Issue, IType, filter_parse
from . import get_orm_params, get_query
from ..services import log_exception, get_bioformat, get_or_create


##
# CREATE
##

def __check_seq_values(**values):
    """
    Process the values 'organism', 'stock', and 'type' into valid IDs,
    and fill in the empty not null fields whenever possible.
    """
    if not values.get('uniquename'):
        raise Exception('Missing the uniquename (sequences ID)')

    if values.get('residues') and not values.get('seqlen'):
        values['seqlen'] = len(values.get('residues'))

    if not values.get('organism_id'):
        org = values.get('organism').strip()
        if org:
            # TODO: check canónical name (getting genus and species?)
            genus = org.split()[0]
            species = org[len(genus):].strip()
            values['organism_id'] = get_or_create(chado_session, Organism, genus=genus, species=species).organism_id
        if not values.get('organism_id'):
            values['organism_id'] = get_or_create(chado_session, Organism, genus='Organism', species='unclassified').organism_id

    if not values.get('type_id') and values.get('type'):
        try:
            values['type_id'] = get_cvterm_query(type=values.get('type'), subtype=values.get('subtype'))[0].one().cvterm_id
        except Exception as e:
            pass
    if not values.get('type_id'):
        try:
            values['type_id'] = get_cvterm_query(type='sequence')[0].one().cvterm_id
        except Exception as e:
            raise Exception(f'Missing the type_id for {values.get("uniquename")}')

    return get_orm_params(Feature, **values)


def __seq_stock(seq, **values):
    if not values.get('stock_id') and values.get('stock'):
        from .individuals import __get_query as get_stock_query, create as create_stock
        try:
            stock = get_stock_query(values.get('stock_id'), uniquename=values.get('stock'))[0].one()
        except Exception as e:
            stock = create_stock(uniquename=values.get('stock'), organism_id=seq.organism_id)[1]
            chado_session.flush()   # stock_id required
        values['stock_id'] = stock.stock_id
    elif not values.get('stock_id'):
        return None
    stock_bind = get_or_create(chado_session, StockFeature,
                               stock_id=values.get('stock_id'),
                               feature_id=seq.feature_id,
                               type_id=seq.type_id)
    return stock


def create(**kwargs):
    content = None
    try:
        values = __check_seq_values(**kwargs)
        content = Feature(**values)
        chado_session.add(content)
        stock = __seq_stock(content, **kwargs)
        # seq to bcs
        bcs_seq = get_or_create(db_session, Sequence,   # specimen_id=bcs_specimen.id
                                chado_id=content.feature_id,
                                chado_table='feature',
                                name=content.uniquename)
        issues, status = [Issue(IType.INFO, f'CREATE sequences: The sequence "{kwargs.get("uniquename")}" was created successfully.')], 201
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, f'CREATE sequences: The sequence "{kwargs.get("uniquename")}" could not be created.')], 409
    return issues, content, status


##
# READ
##

def read(id=None, **kwargs):
    content, count = None, 0
    try:
        content, count = __get_query(id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ sequences: The sequences were successfully read')], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, 'READ sequences: The sequences could not be read.')], 400
    return issues, content, count, status


##
# UPDATE
##

def update(id, **kwargs):
    content = None
    try:
        seq = __get_query(id)[0].one()
        seq.update(get_orm_params(Feature, **kwargs))
        issues, status = [Issue(IType.INFO, f'UPDATE sequences: The sequence "{seq.uniquename}" was successfully updated.')], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, f'UPDATE sequences: The sequence "{seq.uniquename}" could not be updated.')], 409
    return issues, content, status


##
# DELETE
##

def __delete_from_bcs(query):
    ids = [seq.feature_id for seq in query.all()]
    db_session.query(Sequence).filter(Sequence.chado_id.in_(ids)) \
        .delete(synchronize_session='fetch')


def delete(id=None, **kwargs):
    content = None
    try:
        content, count = __get_query(id, **kwargs)
        # TODO: check why is deleting all in bcs with out filter
        __delete_from_bcs(content)
        content = content.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE sequences: {content} sequences were successfully removed.')], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, f'DELETE sequences: The sequences could not be removed.')], 404
    return issues, content, status


##
# IMPORT
##

def __simpleSeq2chado(seq, **params):
    return create(uniquename=seq.id, name=seq.description, residues=str(seq.seq), seqlen=len(seq.seq), **params)[0]


def __fasSeq2chado(seq, **params):
    # params = { organism:<organism>, features=<features>, origin:<origin> }
    features = params.get('features')
    origin = params.get('origin')
    if features:
        if isinstance(features, str):
            features = features.split(',')
        features = [ {'type': f.split()[-1], 'qualifiers': {f.split()[-1]:f[:-len(f.split()[-1])].strip()}} for f in features ]
        params['features'] = features if hasattr(features, '__iter__') else params['features']
        params['molecule_type'] = features[-1]['type']
    elif origin:
        params['molecule_type'] = params.get('origin')
    return create(uniquename=seq.id, name=seq.description, stock=seq.name, residues=str(seq.seq), seqlen=len(seq.seq), **params)[0]


def __gbSeq2chado(seq, **params):
    # seq.id, seq.name, seq.seq, seq.description, seq.dbxrefs, seq.letter_annotations
    # seq.annotations # (<class 'dict'>):
    """
    {
        'molecule_type': 'DNA',
        'topology': 'linear',
        'data_file_division': 'PLN',
        'date': '03-DEC-2018',
        'accessions': ['Seq1_18037'],
        'keywords': [''],
        'source': 'plastid Ruta montana',
        'organism': 'Ruta montana',
        'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'eudicotyledons', 'Gunneridae', 'Pentapetalae', 'rosids', 'malvids', 'Sapindales', 'Rutaceae', 'Rutoideae', 'Ruta'],
        'references': [
            {
                'location': '[0: 1509]',
                'authors': 'Jaen-Molina,R., Soto,M., Marrero,A., Mesa,R. and Caujape-Castells,J.',
                'title': 'Genetic data on the three Canarian endemic Ruta (Rutaceae) provide evidences of overlooked diversification and a complex evolutionary history',
                'journal': 'unpublished',
                'medline_id':'',
                'pubmed_id':'',
                'comment':'',
            }],
        'comment': 'Bankit Comment: ALT EMAIL:rjaenm@grancanaria.com\nBankit Comment: TOTAL # OF SEQS:34',
        'structured_comment': # OrderedDict
            {'Assembly-Data': # OrderedDict
                 {'Sequencing Technology': 'Sanger dideoxy sequencing'}
             }
    }
    """
    # seq.features # (<class 'list'>):
    """
    [
        # .location .type
        SeqFeature(
            location=FeatureLocation(ExactPosition(0), ExactPosition(1509), strand=1),
            qualifiers= # OrderedDict
            {
                'organism': ['Ruta montana'],
                'organelle': ['plastid'],
                'mol_type': ['genomic DNA'],
                'specimen_voucher': ['LPA s.n.'],
                'bio_material': ['JBCVC: DNABANK: 18037'],
                'db_xref': ['taxon:266085'],
                'country': ['Morocco: Between Tanalt to Tidli'],
                'collection_date': ['2015'],
                'collected_by': ['J. Caujape-Castells, C. Harrouni, F. Msanda, et al']
            },
            type='source'),
        SeqFeature(...,
                   qualifiers=OrderedDict({'gene': ['matK']}),
                   type='gene'),
        SeqFeature(...,
                   qualifiers=OrderedDict({'gene': ['matK'],
                                           'note': ['maturase K']}),
                   type='gene'),
        SeqFeature(...,
                   qualifiers=OrderedDict({'gene': ['matK'],
                                           'codon_start': ['1'],
                                           'transl_table': ['11'],
                                           'product': ['maturase K'],
                                           'translation': ['FQVYFELDRSQQHNF...']}),
                   type='CDS')
    ]
    """
    seq.name = f'{seq.id} {seq.description}' if seq.id == seq.name else seq.name
    return create(uniquename=seq.id, name=seq.name, residues=str(seq.seq), seqlen=len(seq.seq),
                  description=seq.description, dbxrefs=seq.dbxrefs, **seq.annotations, features=seq.features, **params)[0]


# ?META: accession, genus, species, molec_region, molec_origin, isle, georegion, individual, collection
def __bio2chado(seq, format, **params):
    try:
        if format == 'fasta':    # fasta
            pattern = '^(?P<id>\w+?) *\[organism=(?P<organism>.+?)\] *(?P<features>.+?); *(?P<origin>.+?)$'
            # p.e.: >Seq1_18037 [organism=Ruta montana] maturase K (matk) gene, partial sequence, partial cds; chloroplast
            meta = re.match(pattern, seq.description)
            if meta:
                return __fasSeq2chado(seq, **meta.groupdict(), **params)
        elif format == 'genbank':   # genbank
            return __gbSeq2chado(seq, **params)
        return __simpleSeq2chado(seq, **params)
    except Exception as e:
        return [Issue(IType.ERROR, f'IMPORT sequences: {seq.id} could not be imported.')]


def import_file(input_file, format=None, **kwargs):
    issues, content, count, status = [], None, 0, 200
    format = get_bioformat(input_file, format)
    try:
        from ..services import seqs_parser
        for s in seqs_parser(input_file, format):
            iss = __bio2chado(s, format, **kwargs)
            issues += [ Issue(i.itype, i.message, os.path.basename(input_file)) for i in iss ]
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR,
                                f'IMPORT sequences: file {input_file} could not be imported.',
                                os.path.basename(input_file))], 409
    return issues, content, status


##
# EXPORT
##

def __seqs_header_parser(seqs, format):     # return dict(uniquename, header)
    # TODO: seqs header parser
    headers = {}
    if format == 'organism' or format == 'organism_canon':
        orgs = chado_session.query(Feature.uniquename, Organism.genus, Organism.species) \
            .join(Organism).filter(Feature.uniquename.in_([x.uniquename for x in seqs])).all()
        from .species_names import get_canonical_species_names
        underscored = True if format == 'organism_canon_underscored' else False
        for seqID, genus, species in orgs:
            headers[seqID] = genus + ' ' + species
            if format in ('organism_canon', 'organism_canon_underscored'):
                headers[seqID] = get_canonical_species_names(db_session,
                                                             [headers[seqID]],
                                                             underscores=underscored)[0] \
                                or headers[seqID].replace(' ', '_') if underscored else headers[seqID]
    return headers


def __seqs2file(seqs, format='fasta', output_file=f"/tmp/output_ngd", header_format=None):
    headers = __seqs_header_parser(seqs, header_format) if header_format else {}
    with open(output_file, "w") as file:
        for seq in seqs:
            if format=='fasta':
                file.write(f'>{headers[seq.uniquename] if headers else seq.uniquename}\n{seq.residues}\n')
            else:
                file.write(f'>{headers[seq.uniquename] if headers else seq.uniquename}\n{seq.residues}\n')
    if format in ('nexus'):
        from Bio import SeqIO
        res=[]
        for s in SeqIO.parse(output_file, 'fasta'):
            s.id = s.description
            s.annotations['molecule_type'] = 'DNA'
            res.append(s)
        SeqIO.write(res, output_file, format=format)
    return output_file


def export(id=None, format='fasta', **kwargs):
    content = None
    try:
        query, count = __get_query(id, **kwargs)
        content = __seqs2file(query.all(), format=format, output_file=f'/tmp/output_ngd.{format}', header_format=kwargs.get('header'))
        issues, status = [Issue(IType.INFO, f'EXPORT sequences: {count} sequences were successfully exported.')], 200
    except Exception as e:
        log_exception(e)
        issues, status = [Issue(IType.ERROR, f'EXPORT sequences: The sequences could not be exported.')], 409
    return issues, content, status


##
# GETTER AND OTHERS
##

def __get_query(id=None, **kwargs):
    if id:
        query = chado_session.query(Feature).filter(Feature.feature_id == id)
        return query, query.count()
    return get_query(chado_session, Feature, aux_filter=__aux_own_filter, aux_order=__aux_own_order, **kwargs)


def __aux_own_filter(filter):
    clause=[]

    if 'analysis_id' in filter:
        from ..db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.feature_id)\
            .filter(filter_parse(AnalysisFeature, [{'analysis_id': filter.get('analysis_id')}]))
        clause.append(Feature.feature_id.in_(_ids))

    if 'phylotree_id' in filter:
        from ..db_models.chado import Phylonode
        _ids = chado_session.query(Phylonode.feature_id)\
            .filter(filter_parse(Phylonode, [{'phylotree_id': filter.get('phylotree_id')}]))
        clause.append(Feature.feature_id.in_(_ids))

    if "prop_cvterm_id" in filter:
        from ..db_models.chado import Featureprop
        _ids = chado_session.query(Featureprop.feature_id)\
            .filter(filter_parse(Featureprop, [{'type_id': filter.get('prop_cvterm_id')}]))
        clause.append(Feature.feature_id.in_(_ids))

    if "program" in filter:
        from ..db_models.chado import Analysis
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, [{'program': filter.get('program')}]))
        from ..db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.feature_id) \
            .filter(AnalysisFeature.analysis_id.in_(_ids))
        clause.append(Feature.feature_id.in_(_ids))

    if "programversion" in filter:
        from ..db_models.chado import Analysis
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, [{'programversion': filter.get('programversion')}]))
        from ..db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.feature_id) \
            .filter(AnalysisFeature.analysis_id.in_(_ids))
        clause.append(Feature.feature_id.in_(_ids))

    if "algorithm" in filter:
        from ..db_models.chado import Analysis
        _ids = chado_session.query(Analysis.analysis_id) \
            .filter(filter_parse(Analysis, [{'algorithm': filter.get('algorithm')}]))
        from ..db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.feature_id) \
            .filter(AnalysisFeature.analysis_id.in_(_ids))
        clause.append(Feature.feature_id.in_(_ids))

    from datetime import datetime
    if "added-from" in filter:
        filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
        _ids = chado_session.query(Feature.feature_id) \
            .filter(filter_parse(Feature, {'timeaccessioned':filter.get("added-from")}))
        clause.append(Feature.feature_id.in_(_ids))
    if "added-to" in filter:
        filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
        _ids = chado_session.query(Feature.feature_id) \
            .filter(filter_parse(Feature, {'timeaccessioned':filter.get("added-to")}))
        clause.append(Feature.feature_id.in_(_ids))
    if "lastmodified-from" in filter:
        filter["lastmodified-from"]['unary'] = datetime.strptime(filter.get("lastmodified-from")['unary'], '%Y-%m-%d')
        _ids = chado_session.query(Feature.feature_id) \
            .filter(filter_parse(Feature, {'timelastmodified':filter.get("lastmodified-from")}))
        clause.append(Feature.feature_id.in_(_ids))
    if "lastmodified-to" in filter:
        filter["lastmodified-to"]['unary'] = datetime.strptime(filter.get("lastmodified-to")['unary'], '%Y-%m-%d')
        _ids = chado_session.query(Feature.feature_id) \
            .filter(filter_parse(Feature, {'timelastmodified':filter.get("lastmodified-to")}))
        clause.append(Feature.feature_id.in_(_ids))

    return clause


def __aux_own_order(order):
    clause=[]
    return clause