from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.db_models.chado import Organism
from biobarcoding.rest import Issue, IType, filter_parse, paginator


def create(**kwargs):
    content = None
    if not kwargs.get('genus'):
        kwargs['genus']='organism'
    if not kwargs.get('species'):
        kwargs['species']='undefined'
    try:
        from biobarcoding.services import get_or_create
        content = get_or_create(chado_session, Organism, **kwargs)
        chado_session.merge(content)
        issues, status = [Issue(IType.INFO, f'CREATE organisms: The organism "{kwargs.get("genus")} {kwargs.get("species")}" was  successfully created.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE organisms: The organism "{kwargs.get("genus")} {kwargs.get("species")}" could not be created.')], 409
    return issues, content, status


count = 0
def read(id = None, **kwargs):
    content = None
    try:
        content = __get_query(id, **kwargs)
        if id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, f'READ organisms: The organisms were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'READ organisms: The organisms could not be read.')], 400
    return issues, content, count, status


def update(id, **kwargs):
    issues = [Issue(IType.WARNING, 'UPDATE phylotrees: dummy completed')]
    return issues, None, 200


def delete(id = None, **kwargs):
    content = None
    try:
        content = __get_query(id, **kwargs).delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, 'DELETE organisms: The organisms were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE organisms: The organisms could not be removed.')], 404
    return issues, content, status


def export(id = None, format = None, output_file = None, **kwargs):
    if not output_file:
        output_file = '/tmp/output_taxa.gbk'
    if not format:
        format = 'genbank'
    try:
        import sys
        stdout = sys.stdout
        with open(output_file, "w") as sys.stdout:
            if format == 'genbank':
                __print_gbk(id, **kwargs)
        sys.stdout = stdout
        issues, status = Issue(IType.INFO, 'EXPORT organisms: The organisms were successfully exported.'), 200
    except Exception as e:
        print(e)
        if stdout:
            sys.stdout = stdout
        issues, status = Issue(IType.ERROR, 'EXPORT organisms: The organisms could not be exported.'), 404
    return issues, output_file, status


def __print_gbk(id=None, **kwargs):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    for org in __get_query(id, **kwargs).all():
        conn.export.export_gbk(org['organism_id'])


def __get_query(id=None, **kwargs):
    query = chado_session.query(Organism)
    global count
    count = 0
    if id:
        query = query.filter(Organism.organism_id==id)
    else:
        if 'filter' in kwargs:
            query = query.filter(filter_parse(Organism, kwargs.get('filter'), __aux_own_filter))
        if 'order' in kwargs:
            query = __get_query_ordered(query, kwargs.get('order'))
        if 'pagination' in kwargs:
            count = query.count()
            query = paginator(query, kwargs.get('pagination'))
    return query


def __aux_own_filter(filter):
    clause = []
    if filter.get('feature_id'):
        from biobarcoding.db_models.chado import Feature
        _organism_ids = chado_session.query(Feature.organism_id)\
            .filter(filter_parse(Feature, [{'feature_id': filter.get('feature_id')}])).all()
        clause.append(Organism.organism_id.in_(_organism_ids))
    if filter.get('analysis_id'):
        from biobarcoding.db_models.chado import AnalysisFeature, Feature
        _feature_ids = chado_session.query(AnalysisFeature.feature_id)\
            .filter(filter_parse(AnalysisFeature, [{'analysis_id': filter.get('analysis_id')}])).all()
        _organism_ids = chado_session.query(Feature.organism_id)\
            .filter(Feature.feature_id.in_(_feature_ids)).all()
        clause.append(Organism.organism_id.in_(_organism_ids))
    if filter.get('phylonode_id'):
        from biobarcoding.db_models.chado import PhylonodeOrganism
        _organism_ids = chado_session.query(PhylonodeOrganism.organism_id)\
            .filter(PhylonodeOrganism.phylonode_id == filter.get('phylonode_id')).all()
        clause.append(Organism.organism_id.in_(_organism_ids))
    if filter.get('phylotree_id'):
        from biobarcoding.db_models.chado import Phylonode, PhylonodeOrganism
        _phylonode_ids = chado_session.query(Phylonode.phylonode_id)\
            .filter(filter_parse(Phylonode, [{'phylotree_id': filter.get('phylotree_id')}])).all()
        _organism_ids = chado_session.query(PhylonodeOrganism.organism_id)\
            .filter(PhylonodeOrganism.phylonode_id.in_(_phylonode_ids)).all()
        clause.append(Organism.organism_id.in_(_organism_ids))
    if filter.get('rank'):
        # TODO: organism.type_id ?
        from biobarcoding.db_models.chado import Phylonode, PhylonodeOrganism, Cv, Cvterm
        _rank_cvterm_ids = chado_session.query(Cvterm.cvterm_id)\
            .join(Cv).filter(Cv.name=='taxonomy')\
            .filter(filter_parse(Cvterm, [{'name':filter.get('rank')}])).all()
        _phylonode_ids = chado_session.query(Phylonode.phylonode_id)\
            .filter(Phylonode.type_id.in_(_rank_cvterm_ids)).all()
        _ids = chado_session.query(PhylonodeOrganism.organism_id)\
            .filter(PhylonodeOrganism.phylonode_id.in_(_phylonode_ids)).all()
        clause.append(Organism.organism_id.in_(_ids))
    return clause


def __get_query_ordered(query, order):
    # query = query.order(order_parse(Organism, kwargs.get('order'), __aux_own_order))
    return query
