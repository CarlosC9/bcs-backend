from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.rest import Issue, IType


def create_organisms(genus, species, common_name = None, abbreviation = None, comment = None):
    content = { 'genus':genus, 'species':species, 'common_name':common_name, 'abbreviation':abbreviation, 'comment':comment }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        from biobarcoding.services import get_or_create
        from biobarcoding.db_models.chado import Organism
        content = get_or_create(chado_session, Organism,
            genus=genus,
            species=species,
            common_name=common_name,
            abbreviation=abbreviation,
            comment=comment)
        chado_session.merge(content)
        issues, status = [Issue(IType.INFO, f'CREATE organisms: The organism "{genus} {species}" was  successfully created.')], 201
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'CREATE organisms: The organism "{genus} {species}" could not be created.')], 500
    return issues, content, status


def read_organisms(organism_id = None, genus = None, species = None, common_name = None, abbreviation = None, comment = None):
    content = { 'organism_id':organism_id, 'genus':genus, 'species':species, 'common_name':common_name, 'abbreviation':abbreviation, 'comment':comment }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        content = __get_query(organism_id, genus=genus, species=species,
                common_name=common_name, abbreviation=abbreviation, comment=comment)
        if organism_id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, f'READ organisms: The organisms were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'READ organisms: The organisms could not be read.')], 500
    return issues, content, status


def update_organisms(organism_id, genus = None, species = None, common_name = None, abbreviation = None, comment = None):
    issues = [Issue(IType.WARNING, 'UPDATE phylotrees: dummy completed')]
    content = { 'organism_id':organism_id, 'genus':genus, 'species':species, 'common_name':common_name, 'abbreviation':abbreviation, 'comment':comment }
    return issues, content, 200


def delete_organisms(organism_id = None, ids = None, genus = None, species = None, common_name = None, abbreviation = None, comment = None):
    content = { 'organism_id':organism_id, 'genus':genus, 'species':species, 'common_name':common_name, 'abbreviation':abbreviation, 'comment':comment }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        content = __get_query(organism_id, ids, genus, species,
            common_name, abbreviation, comment).delete(synchronize_session='fetch')
        issues, status = [Issue(IType.INFO, 'DELETE organisms: The organisms were successfully removed.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'DELETE organisms: The organisms could not be removed.')], 500
    return issues, content, status


def export_organisms(organism_id = None, format = 'genbank', output_file = None):
    if not output_file:
        output_file = '/tmp/output_taxa.gbk'
    try:
        from biobarcoding.services import conn_chado
        conn = conn_chado()
        import sys
        stdout = sys.stdout
        with open(output_file, "w") as sys.stdout:
            if organism_id:
                conn.export.export_gbk(organism_id)
            else:
                orgs = conn.organism.get_organisms()
                for org in orgs:
                    conn.export.export_gbk(org['organism_id'])
        sys.stdout = stdout
        issues, status = Issue(IType.INFO, 'EXPORT organisms: The organisms were successfully exported.'), 200
    except Exception as e:
        print(e)
        if stdout:
            sys.stdout = stdout
        issues, status = Issue(IType.ERROR, 'EXPORT organisms: The organisms could not be exported.'), 500
    return issues, output_file, status


def __get_query(organism_id=None, ids=None, genus=None, species=None, common_name=None, abbreviation=None, comment=None, feature_id=None, phylonode_id=None):
    from biobarcoding.db_models.chado import Organism
    query = chado_session.query(Organism)
    if organism_id:
        query = query.filter(Organism.organism_id==organism_id)
    if ids:
        query = query.filter(Organism.organism_id.in_(ids))
    if genus:
        query = query.filter(Organism.genus==genus)
    if species:
        query = query.filter(Organism.species==species)
    if abbreviation:
        query = query.filter(Organism.abbreviation==abbreviation)
    if common_name:
        query = query.filter(Organism.common_name==common_name)
    if comment:
        query = query.filter(Organism.comment==comment)
    if feature_id:
        from biobarcoding.db_models.chado import Feature
        feat_org_id = chado_session.query(Feature.organism_id)\
            .filter(Feature.feature_id==feature_id)
        query = query.filter(Organism.organism_id==feat_org_id)
    if phylonode_id:
        from biobarcoding.db_models.chado import PhylonodeOrganism
        node_org_id = chado_session.query(PhylonodeOrganism.organism_id)\
            .filter(PhylonodeOrganism.phylonode_id==phylonode_id)
        query = query.filter(Organism.organism_id==node_org_id)
    return query
