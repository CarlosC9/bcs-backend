from biobarcoding.db_models import DBSession as db_session
from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.db_models.chado import Feature

from biobarcoding.rest import Issue, IType, filter_parse, paginator


def create(**kwargs):
    issues = [Issue(IType.WARNING, 'CREATE sequences: dummy completed')]
    return issues, None, 200


def read(id=None, **kwargs):
    content = None
    try:
        content = __get_query(id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ sequences: The sequences were successfully read')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ sequences: The sequences could not be read.')], 500
    return issues, content, status


def update(id, **kwargs):
    issues = [Issue(IType.WARNING, 'UPDATE sequences: dummy completed')]
    return issues, None, 200


def __delete_from_bcs(feature_id):
    from biobarcoding.db_models.bioinformatics import Sequence
    db_session.query(Sequence).filter(Sequence.chado_feature_id == feature_id) \
        .delete(synchronize_session='fetch')


def delete(id=None, **kwargs):
    content = None
    try:
        query = __get_query(id, **kwargs)
        for seq in query.all():
            __delete_from_bcs(seq.feature_id)
        resp = query.delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, f'DELETE sequences: {resp} sequences were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'DELETE sequences: The sequences could not be removed.')], 500
    return issues, content, status

# TODO: replace python-chado lib?
def import_file(input_file, format='fasta', **kwargs):
    content = None
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if not kwargs.get('organism_id'):
        try:
            kwargs['organism_id'] = conn.organism.add_organism(genus='organism',
                                                     species='undefined', common='', abbr='')['organism_id']
        except Exception as e:
            kwargs['organism_id'] = conn.organism.get_organisms(species='undefined')[0]['organism_id']
    try:
        resp = conn.feature.load_fasta(input_file, kwargs.get('organism_id'), analysis_id=kwargs.get('analysis_id'), update=True)
        from Bio import SeqIO
        __seqs2bcs([seq.id for seq in SeqIO.parse(input_file, format or 'fasta')])
        issues, status = [Issue(IType.INFO, f'IMPORT sequences: {resp} sequences were successfully imported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'IMPORT sequences: file {input_file} could not be imported.')], 500
    return issues, content, status


def __seqs2bcs(names):
    from biobarcoding.db_models.chado import Feature
    from biobarcoding.db_models import DBSession as db_session
    for seq in chado_session.query(Feature).filter(Feature.uniquename.in_(names)).all():
        db_session.merge(__feature2bcs(seq))


def __feature2bcs(seq):
    from biobarcoding.services import get_or_create
    from biobarcoding.db_models.bioinformatics import Specimen, Sequence
    bcs_specimen = get_or_create(db_session, Specimen, name=seq.uniquename)
    db_session.merge(bcs_specimen)
    db_session.flush()
    bcs_sequence = get_or_create(db_session, Sequence,
                                 chado_feature_id=seq.feature_id,
                                 chado_table='feature',
                                 name=seq.uniquename,
                                 specimen_id=bcs_specimen.id)
    return bcs_sequence


def __seqs2file(output_file, format, seqs):
    with open(output_file, "w") as file:
        for seq in seqs:
            if format=='fasta':
                file.write(f'>{seq.uniquename}\n{seq.residues}\n')
            else:
                file.write(f'>{seq.uniquename}\n{seq.residues}\n')


def export(id=None, output_file="/tmp/output_seqs.fas", format='fasta', **kwargs):
    try:
        query = __get_query(id, **kwargs)
        __seqs2file(output_file, format, query.all())
        issues, status = [Issue(IType.INFO, f'EXPORT sequences: {query.count()} sequences were successfully exported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'EXPORT sequences: The sequences could not be exported.')], 500
    return issues, output_file, status


def __get_query(id=None, **kwargs):
    query = chado_session.query(Feature)
    if id:
        query = query.filter(Feature.feature_id == id)
    else:
        if 'filter' in kwargs:
            query = query.filter(filter_parse(Feature, kwargs.get('filter'), __aux_own_filter))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            query = paginator(query, kwargs.get('pagination'))
    return query


def __aux_own_filter(filter):
    clause=[]
    if 'analysis_id' in filter:
        from biobarcoding.db_models.chado import AnalysisFeature
        _ids = chado_session.query(AnalysisFeature.feature_id)\
                .filter(
                    filter_parse(AnalysisFeature, [{'analysis_id': filter.get('analysis_id')}]))
        clause.append(Feature.feature_id.in_(_ids))
    if 'phylotree_id' in filter:
        from biobarcoding.db_models.chado import Phylonode
        _ids = chado_session.query(Phylonode.feature_id)\
                .filter(
                    filter_parse(Phylonode, [{'phylotree_id': filter.get('phylotree_id')}]))
        clause.append(Feature.feature_id.in_(_ids))
    return clause


def __get_query_ordered(query, order):
    # query = query.order(order_parse(Feature, kwargs.get('order'), __aux_own_order))
    return query