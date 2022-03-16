from ....db_models import DBSession as db_session, DBSessionChado as chado_session
from ....db_models.chado import Organism
from ....rest import Issue, IType, filter_parse
from ....common import generate_json
from ... import get_query, force_underscored


def create(**kwargs):
    content = None
    if not kwargs.get('genus'):
        kwargs['genus'] = 'Organism'
    if not kwargs.get('species'):
        kwargs['species'] = 'Unclassified'
    try:
        from ... import get_or_create
        content = get_or_create(chado_session, Organism, **kwargs)
        issues, status = [Issue(IType.INFO, f'CREATE organisms: The organism "{kwargs.get("genus")} {kwargs.get("species")}" was  successfully created.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE organisms: The organism "{kwargs.get("genus")} {kwargs.get("species")}" could not be created.')], 409
    return issues, content, status


def get_canonical(organism_id=None, name=None, underscores=False, **kwargs):
    if organism_id and not name:
        name = chado_session.query(Organism.name)\
            .filter(Organism.organism_id == organism_id).one()
    from ...species_names import get_canonical_species_names
    return get_canonical_species_names(db_session,
                                       [name],
                                       underscores=underscores)[0] \
           or force_underscored(name) if underscores else name


def __append_canonical(*orgs):
    res = []
    for org in orgs:
        import json
        new = json.loads(generate_json(org))
        new['name'] = "%s %s" % (new.get('genus'), new.get('species'))
        new['canonical_name'] = get_canonical(**new)
        new['canonical_underscored_name'] = get_canonical(**new, underscores=True)
        res.append(new)
    return res


def read(id=None, **kwargs):
    content, count = None, 0
    try:
        content, count = __get_query(id, **kwargs)
        if id:
            content = __append_canonical(content.one())
        else:
            content = __append_canonical(*content.all())
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
        content, count = __get_query(id, **kwargs)
        content = content.delete(synchronize_session='fetch')
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
    from ... import conn_chado
    conn = conn_chado()
    for org in __get_query(id, **kwargs)[0].all():
        conn.export.export_gbk(org['organism_id'])


def __get_query(organism_id=None, **kwargs):
    query = chado_session.query(Organism)
    if organism_id:
        query = query.filter(Organism.organism_id == organism_id)
        return query, query.count()
    return get_query(chado_session, Organism, query=query, aux_filter=__aux_own_filter, aux_order=__aux_own_order, **kwargs)


def __aux_own_filter(filter):
    clause = []
    if filter.get('feature_id'):
        from ....db_models.chado import Feature
        _organism_ids = chado_session.query(Feature.organism_id)\
            .filter(filter_parse(Feature, [{'feature_id': filter.get('feature_id')}])).all()
        clause.append(Organism.organism_id.in_(_organism_ids))
    if filter.get('analysis_id'):
        from ....db_models.chado import AnalysisFeature, Feature
        _feature_ids = chado_session.query(AnalysisFeature.feature_id)\
            .filter(filter_parse(AnalysisFeature, [{'analysis_id': filter.get('analysis_id')}])).all()
        _organism_ids = chado_session.query(Feature.organism_id)\
            .filter(Feature.feature_id.in_(_feature_ids)).all()
        clause.append(Organism.organism_id.in_(_organism_ids))
    if filter.get('phylonode_id'):
        from ....db_models.chado import PhylonodeOrganism
        _organism_ids = chado_session.query(PhylonodeOrganism.organism_id)\
            .filter(PhylonodeOrganism.phylonode_id == filter.get('phylonode_id')).all()
        clause.append(Organism.organism_id.in_(_organism_ids))
    if filter.get('phylotree_id'):
        from ....db_models.chado import Phylonode, PhylonodeOrganism
        _phylonode_ids = chado_session.query(Phylonode.phylonode_id)\
            .filter(filter_parse(Phylonode, [{'phylotree_id': filter.get('phylotree_id')}])).all()
        _organism_ids = chado_session.query(PhylonodeOrganism.organism_id)\
            .filter(PhylonodeOrganism.phylonode_id.in_(_phylonode_ids)).all()
        clause.append(Organism.organism_id.in_(_organism_ids))
    if filter.get('rank'):
        # TODO: organism.type_id ?
        from ....db_models.chado import Phylonode, PhylonodeOrganism, Cv, Cvterm
        _rank_cvterm_ids = chado_session.query(Cvterm.cvterm_id)\
            .join(Cv).filter(Cv.name=='taxonomy')\
            .filter(filter_parse(Cvterm, [{'name':filter.get('rank')}])).all()
        _phylonode_ids = chado_session.query(Phylonode.phylonode_id)\
            .filter(Phylonode.type_id.in_(_rank_cvterm_ids)).all()
        _ids = chado_session.query(PhylonodeOrganism.organism_id)\
            .filter(PhylonodeOrganism.phylonode_id.in_(_phylonode_ids)).all()
        clause.append(Organism.organism_id.in_(_ids))
    return clause


def __aux_own_order(order):
    # query = query.order(order_parse(Organism, kwargs.get('order'), __aux_own_order))
    return []
