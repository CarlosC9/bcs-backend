def create_organisms(genus, species, common = None, abbr = None, comment = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    try:
        resp = conn.organism.add_organism(
            genus=genus,
            species=species,
            common=common,
            abbr=abbr,
            comment=comment)
        return {'status':'success','message':f'The organism "{genus} {species}" created successfully.'}, 201
    except Exception as e:
        return {'status':'failure','message':f'The organism "{genus} {species}" could not be created.'}, 500


def read_organisms(organism_id = None, genus = None, species = None, common = None, abbr = None, comment = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    resp = conn.organism.get_organisms(organism_id, genus, species, common, abbr, comment)
    if organism_id:
        return resp[0], 200
    return resp, 200


def update_organisms(organism_id, genus = None, species = None, common = None, abbr = None, comment = None):
    return {'status':'success','message':'UPDATE: organisms dummy completed'}, 200


def delete_organisms(organism_id = None, ids = None, genus = None, species = None, common = None, abbr = None, comment = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    resp = 0
    try:
        if not ids:
            resp = conn.organism.delete_organisms(organism_id, genus, species, common, abbr, comment)
        else:
            for id in ids:
                resp += conn.organism.delete_organisms(id, genus, species, common, abbr, comment)
        return {'status':'success','message':f'{resp} organisms were successfully removed.'}, 200
    except Exception as e:
        return {'status':'success','message':f'{resp} organisms were removed.'}, 207
        # return {'status':'failure','message':f'The organisms could not be removed.'}, 500


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
