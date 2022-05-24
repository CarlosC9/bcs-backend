from typing import List

from . import MetaService
from ...main import get_orm
from ....db_models import DBSessionChado
from ... import log_exception, get_or_create

NAMES = {
    'gbif': 'GBIF taxonomy tree',
    'ncbi': 'NCBI taxonomy tree',
    'biota': 'BIOTA taxonomy tree',
}

##
# TAXONOMY SERVICE
##
class Service(MetaService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('taxonomies')

    ##
    # CREATE
    ##

    def check_values(self, **values) -> dict:
        if not values.get('dbxref_id'):
            dbxref = values.get('dbxref', f'taxonomy:{values.get("name")}')
            if not dbxref.startswith('taxonomy:'):
                dbxref = 'taxonomy:' + dbxref
            import time
            from ....db_models.chado import Db, Dbxref
            values['dbxref_id'] = get_or_create(self.db, Dbxref,
                                                db_id=self.db.query(Db).filter(Db.name == 'null').one().db_id,
                                                accession=dbxref,
                                                version=time.strftime("%Y %b %d %H:%M:%S")).dbxref_id
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

    def get_query(self, **kwargs):
        query = None
        from ....db_models.chado import Dbxref
        query = self.db.query(self.orm).filter(self.orm.dbxref_id.in_(
            self.db.query(Dbxref.dbxref_id).filter(Dbxref.accession.like('taxonomy:%')).subquery()))
        return super(Service, self).get_query(query=query, **kwargs)

    def aux_filter(self, filter):
        from ....rest import filter_parse
        clauses = []

        if filter.get('organism_id'):
            from ....db_models.chado import Phylonode, PhylonodeOrganism
            _ids = self.db.query(Phylonode.phylotree_id).filter(Phylonode.phylonode_id.in_(
                self.db.query(PhylonodeOrganism.phylonode_id).filter(
                    filter_parse(PhylonodeOrganism, [{'organism_id': filter.get('organism_id')}])).subquery())
            ).subquery()
            clauses.append(self.orm.phylotree_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(filter)


def insert_taxon(**kwargs) -> List[str]:
    return []
