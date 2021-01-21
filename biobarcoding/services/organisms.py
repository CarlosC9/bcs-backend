from biobarcoding.authentication import bcs_session
from biobarcoding.db_models import DBSession as db_session
from biobarcoding.db_models import DBSessionChado as chado_session


@bcs_session(read_only=False)
def create_organisms(genus, species, common_name = None, abbreviation = None, comment = None):
    from biobarcoding.services import get_or_create
    from biobarcoding.db_models.chado import Organism
    try:
        resp = get_or_create(chado_session, Organism,
            genus=genus,
            species=species,
            common_name=common_name,
            abbreviation=abbreviation,
            comment=comment)
        chado_session.merge(resp)
        return {'status':'success','message':f'The organism "{genus} {species}" created successfully.'}, 201
    except Exception as e:
        return {'status':'failure','message':f'The organism "{genus} {species}" could not be created.'}, 500


@bcs_session(read_only=True)
def read_organisms(organism_id = None, genus = None, species = None, common_name = None, abbreviation = None, comment = None):
    resp = __get_query(organism_id, genus=genus, species=species,
            common_name=common_name, abbreviation=abbreviation, comment=comment)
    from biobarcoding.services import chado2json
    resp = chado2json(resp)
    if organism_id:
        return resp[0], 200
    return resp, 200


@bcs_session(read_only=False)
def update_organisms(organism_id, genus = None, species = None, common_name = None, abbreviation = None, comment = None):
    try:
        resp = __get_query(organism_id, genus=genus, species=species,
            common_name=common_name, abbreviation=abbreviation, comment=comment)
        chado_session.merge(resp)
        return {'status':'success','message':f'The organism "{genus} {species}" updated successfully.'}, 201
    except Exception as e:
        return {'status':'failure','message':f'The organism "{genus} {species}" could not be updated.'}, 500


@bcs_session(read_only=False)
def delete_organisms(organism_id = None, ids = None, genus = None, species = None, common_name = None, abbreviation = None, comment = None):
    try:
        resp = __get_query(organism_id, ids, genus, species,
            common_name, abbreviation, comment).delete(synchronize_session='fetch')
        return {'status':'success','message':f'{resp} organisms were successfully removed.'}, 200
    except Exception as e:
        return {'status':'success','message':f'{resp} organisms were removed.'}, 207
        # return {'status':'failure','message':f'The organisms could not be removed.'}, 500


@bcs_session(read_only=True)
def export_organisms(organism_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    import sys
    stdout = sys.stdout
    with open('/tmp/output_taxa.gbk', "w") as sys.stdout:
        if organism_id:
            conn.export.export_gbk(organism_id)
        else:
            orgs = conn.organism.get_organisms()
            for org in orgs:
                try:
                    conn.export.export_gbk(org['organism_id'])
                except Exception as e:
                    pass
    sys.stdout = stdout
    return '/tmp/' + output_file, 200


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
        query = query.join(Feature)\
            .filter(Feature.feature_id==feature_id)
    if phylonode_id:
        from biobarcoding.db_models.chado import PhylonodeOrganism
        query = query.join(PhylonodeOrganism) \
            .filter(PhylonodeOrganism.phylonode_id==phylonode_id)
    return query
