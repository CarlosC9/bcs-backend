import os.path

from . import MetaService
from ...main import get_orm
from ....db_models import DBSession, DBSessionChado
from ....rest import filter_parse
from ....common import generate_json
from ... import log_exception


##
# ORGANISM SERVICE
##
class Service(MetaService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('organisms')

    ##
    # CREATE
    ##

    def prepare_values(self, **values) -> dict:
        values['species'] = values.get('species', values.get('organism'))
        return super(Service, self).prepare_values(**values)

    def check_values(self, **values) -> dict:
        values['genus'] = values.get('genus', 'unknown')
        values['species'] = values.get('species', 'Unclassified organism')
        return super(Service, self).check_values(**values)

    ##
    # READ
    ##

    def attach_data(self, content):
        import json
        new = json.loads(generate_json(content))
        new['name'] = self.db.query(self.orm.name).filter(self.orm.organism_id == content.organism_id).one()[0]
        from ... import force_underscored
        from ...species_names import get_canonical_species_names
        new['canonical_name'] = get_canonical_species_names(DBSession, new.get('name'))[0] or new.get('name')
        new['canonical_underscored_name'] = get_canonical_species_names(DBSession, new.get('name'), underscores=True)[0] \
                                            or force_underscored(new.get('name'))
        return new

    ##
    # EXPORT
    ##

    def export_file(self, format=None, output_file=None, **kwargs):
        count = 0
        output_file = output_file or '/tmp/output_taxa.gbk'
        format = format or 'genbank'
        try:
            import sys
            stdout = sys.stdout
            with open(output_file, "w") as sys.stdout:
                if format == 'genbank':
                    content, count = self.print_gbk(**kwargs)
            sys.stdout = stdout
        except Exception as e:
            if stdout:
                sys.stdout = stdout
            # issues, status = Issue(IType.ERROR, 'EXPORT organisms: The organisms could not be exported.'), 404
            log_exception(e)
            raise Exception(f'EXPORT organism: file {os.path.basename(output_file)} could not be exported.')
        return output_file, count

    def print_gbk(self, id=None, **kwargs):
        from ... import conn_chado
        conn = conn_chado()
        content, count = self.get_query(id, **kwargs)
        for org in content.all():
            conn.export.export_gbk(org['organism_id'])
        return content, count

    ##
    # GET SQLALCHEMY QUERY
    ##

    def aux_filter(self, filter):
        clauses = []

        if filter.get('feature_id'):
            from ....db_models.chado import Feature
            _organism_ids = self.db.query(Feature.organism_id)\
                .filter(filter_parse(Feature, [{'feature_id': filter.get('feature_id')}])).all()
            clauses.append(self.orm.organism_id.in_(_organism_ids))

        if filter.get('analysis_id'):
            from ....db_models.chado import AnalysisFeature, Feature
            _feature_ids = self.db.query(AnalysisFeature.feature_id)\
                .filter(filter_parse(AnalysisFeature, [{'analysis_id': filter.get('analysis_id')}])).all()
            _organism_ids = self.db.query(Feature.organism_id)\
                .filter(Feature.feature_id.in_(_feature_ids)).all()
            clauses.append(self.orm.organism_id.in_(_organism_ids))

        if filter.get('phylonode_id'):
            from ....db_models.chado import PhylonodeOrganism
            _organism_ids = self.db.query(PhylonodeOrganism.organism_id)\
                .filter(PhylonodeOrganism.phylonode_id == filter.get('phylonode_id')).all()
            clauses.append(self.orm.organism_id.in_(_organism_ids))

        if filter.get('phylotree_id'):
            from ....db_models.chado import Phylonode, PhylonodeOrganism
            _phylonode_ids = self.db.query(Phylonode.phylonode_id)\
                .filter(filter_parse(Phylonode, [{'phylotree_id': filter.get('phylotree_id')}])).all()
            _organism_ids = self.db.query(PhylonodeOrganism.organism_id)\
                .filter(PhylonodeOrganism.phylonode_id.in_(_phylonode_ids)).all()
            clauses.append(self.orm.organism_id.in_(_organism_ids))

        if filter.get('rank'):
            # TODO: organism.type_id ?
            from ....db_models.chado import Phylonode, PhylonodeOrganism, Cv, Cvterm
            _rank_cvterm_ids = self.db.query(Cvterm.cvterm_id)\
                .join(Cv).filter(Cv.name=='taxonomy')\
                .filter(filter_parse(Cvterm, [{'name':filter.get('rank')}])).all()
            _phylonode_ids = self.db.query(Phylonode.phylonode_id)\
                .filter(Phylonode.type_id.in_(_rank_cvterm_ids)).all()
            _ids = self.db.query(PhylonodeOrganism.organism_id)\
                .filter(PhylonodeOrganism.phylonode_id.in_(_phylonode_ids)).all()
            clauses.append(self.orm.organism_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(filter)
