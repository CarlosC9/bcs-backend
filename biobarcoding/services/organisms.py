def create_organisms(genus, species, common = None, abbr = None, comment = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    resp = conn.organism.add_organism(
        genus=genus,
        species=species,
        common=common,
        abbr=abbr,
        comment=comment)
    return {'status':'success','message':resp}, 200


def read_organisms(organism_id = None, genus = None, species = None, common = None, abbr = None, comment = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    resp = conn.organism.get_organisms(organism_id, genus, species, common, abbr, comment)
    return {'status':'success','message':resp}, 200


def update_organisms(organism_id, genus = None, species = None, common = None, abbr = None, comment = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    # resp = conn.organism.update_organism(
    #     genus=genus,
    #     species=species,
    #     common=common,
    #     abbr=abbr)
    return {'status':'success','message':resp}, 200


def delete_organisms(organism_id = None, genus = None, species = None, common = None, abbr = None, comment = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    resp = conn.organism.delete_organisms(organism_id, genus, species, common, abbr, comment)
    return {'status':'success','message':resp}, 200


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
