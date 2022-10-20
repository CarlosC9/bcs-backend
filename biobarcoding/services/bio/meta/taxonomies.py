from typing import List

from sqlalchemy import or_

from . import MetaService
from ...main import get_orm
from ... import log_exception, get_or_create
from ....db_models.chado import Phylonode, PhylonodeOrganism, Organism


##
# TAXONOMY SERVICE
##
class Service(MetaService):

    def __init__(self):
        super(Service, self).__init__()
        self.orm = get_orm('taxonomies')

    ##
    # CREATE
    ##

    def check_values(self, **values) -> dict:

        if not values.get('dbxref_id'):
            dbxref = values.get('dbxref', f'taxonomy:{values.get("name")}')
            if not dbxref.startswith('taxonomy:'):
                dbxref = 'taxonomy:' + dbxref
            # import time
            from ....db_models.chado import Db, Dbxref
            values['dbxref_id'] = get_or_create(self.db, Dbxref,
                                                db_id=self.db.query(Db).filter(Db.name == 'null').one().db_id,
                                                # version=time.strftime("%Y %b %d %H:%M:%S")
                                                accession=dbxref).dbxref_id

        if not values.get('type_id'):
            from .ontologies import get_type_id
            try:
                values['type_id'] = get_type_id(type=values.get('type', 'taxonomy'))
            except:
                pass

        return super(Service, self).check_values(**values)

    ##
    # IMPORT
    ##
    def import_file(self, infile, **kwargs):
        # format: ncbi, gbif, biota ?
        # source: biota ?
        content, count = None, 0
        try:
            from flask import current_app
            cfg = current_app.config
            named = f' -n {kwargs.get("name")} ' if kwargs.get("name") else ''
            import os.path
            dir_path = os.path.dirname(os.path.realpath(__file__))
            from ... import exec_cmds
            content, err = exec_cmds(f'''(cd {dir_path}/perl_scripts/ &&
                perl ./load_ncbi_taxonomy.pl\
                    -H {cfg["CHADO_HOST"]}\
                    -D {cfg["CHADO_DATABASE"]}\
                    -u {cfg["CHADO_USER"]}\
                    -p {cfg["CHADO_PASSWORD"]}\
                    -d Pg\
                    -i {infile}\
                    {named})''')
            self.db.execute(
                "SELECT setval('phylonode_phylonode_id_seq', (SELECT MAX(phylonode_id) FROM phylonode)+1);")
            # self.db.execute("ALTER SEQUENCE phylonode_phylonode_id_seq RESTART WITH (SELECT MAX(phylonode_id) FROM phylonode)+1;")
            if err:
                e = Exception(f'IMPORT taxonomies: The taxonomy {os.path.basename(infile)} was barely imported.')
                log_exception(e)
                # raise e
        except Exception as e:
            log_exception(e)
            raise Exception(f'IMPORT taxonomies: The taxonomy {os.path.basename(infile)} could not be imported.')
        return content, count

    ##
    # GET SQLALCHEMY QUERY
    ##

    def get_query(self, query=None, **kwargs):
        query = query or self.db.query(self.orm)
        from ....db_models.chado import Dbxref      # TODO filter by type_id instead
        query = query.filter(self.orm.dbxref_id.in_(
            self.db.query(Dbxref.dbxref_id).filter(Dbxref.accession.like('taxonomy:%')).subquery()))
        return super(Service, self).get_query(query=query, **kwargs)

    def aux_filter(self, _filter: dict):
        from ....rest import filter_parse
        clauses = []

        if _filter.get('organism_id'):
            _ids = self.db.query(Phylonode.phylotree_id).join(PhylonodeOrganism).filter(
                    filter_parse(PhylonodeOrganism, [{'organism_id': _filter.get('organism_id')}]))
            clauses.append(self.orm.phylotree_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(_filter)


def insert_taxon(organism_id, **kwargs) -> List[str]:

    if kwargs.get('taxonomy_id') or kwargs.get('taxonomy'):
        self = Service()

        def get_phylonode(q_filter):
            return self.db.query(Phylonode).join(self.orm).join(PhylonodeOrganism).join(Organism) \
                .filter(self.orm.phylotree_id == taxonomy.phylotree_id, q_filter).one()

        taxonomy = get_or_create(self.db, self.orm, **self.check_values(
            id=kwargs.get('taxonomy_id'), name=kwargs.get('taxonomy')))

        try:
            node = get_phylonode(PhylonodeOrganism.organism_id == organism_id)
        except Exception as e:
            node = Phylonode(phylotree_id=taxonomy.phylotree_id)
            self.db.add(node)
            self.db.flush()
            get_or_create(self.db, PhylonodeOrganism, phylonode_id=node.phylonode_id, organism_id=organism_id)

        if kwargs.get('parent_id') or kwargs.get('parent'):
            parent = get_phylonode(or_(PhylonodeOrganism.organism_id == kwargs.get('parent_id'),
                                       Organism.species == kwargs.get('parent')))
            node.parent_phylonode_id = parent.phylonode_id

    return []
