def create_organisms(genus, species, common = None, abbr = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    res = conn.organism.add_organism(
        genus=genus,
        species=species,
        common=common,
        abbr=abbr)
    return res


def read_organisms(organism_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if organism_id:
        res = conn.organism.get_organisms(organism_id = organism_id)
    else:
        res = conn.organism.get_organisms()
    return res


def update_organisms(organism_id, genus, species, common = None, abbr = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    # res = conn.organism.update_organism(
    #     genus=genus,
    #     species=species,
    #     common=common,
    #     abbr=abbr)
    return res


def delete_organisms(organism_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if organism_id:
        res = conn.organism.delete_organisms(organism_id = organism_id)
    else:
        res = conn.organism.delete_organisms()
    return res


def export_organisms(output_file = None, organism_id = None):
    from biobarcoding.services import conn_chado
    conn = conn_chado()
    if not output_file:
        output_file = 'output_taxa.gbk'
    import sys
    with open('/tmp/' + output_file, "w") as sys.stdout:
        if organism_id:
            conn.export.export_gbk(organism_id)
        else:
            orgs = conn.organism.get_organisms()
            for org in orgs:
                try:
                    conn.export.export_gbk(org['organism_id'])
                except Exception as e:
                    pass
    return '/tmp/' + output_file
