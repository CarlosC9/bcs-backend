from . import MetaService
from ...main import get_orm
from ....db_models import DBSession, DBSessionChado
from ....rest import filter_parse


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
        new = super(Service, self).attach_data(content)

        if new:
            try:
                new['name'] = " ".join([new['genus'], new['species'], new['infraspecific_name']]).strip()
                from ... import force_underscored
                from ...species_names import get_canonical_species_names
                new['canonical_name'] = get_canonical_species_names(DBSession, [new.get('name')])[0] or new.get('name')
                new['canonical_underscored_name'] = get_canonical_species_names(DBSession, [new.get('name')], underscores=True)[0] \
                                                    or force_underscored(new.get('name'))
            except:
                pass

        return new

    ##
    # EXPORT
    ##

    def data2file(self, orgs: list, outfile, format: str, **kwargs) -> int:
        from Bio import SeqIO
        return SeqIO.write(self.chado2biopy(orgs), outfile, format)

    def chado2biopy(self, orgs: list) -> list:
        from ..meta.ontologies import CvtermService
        from ....db_models.chado import Feature, Featureloc, Featureprop
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        seqs = []
        for org in orgs:
            seq = Seq('')
            record_features = []
            features = self.db.query(Feature, Featureloc) \
                .filter_by(organism_id=org.organism_id) \
                .join(Featureloc, Feature.feature_id == Featureloc.feature_id, isouter=True)

            for idx, (feature, featureloc) in enumerate(features):
                # Sequence containing feature
                if feature.residues:
                    ## This seems bad? What if multiple things have seqs?
                    # seq = Seq(feature.residues)
                    pass
                else:
                    qualifiers = {
                        CvtermService.get_query(purpose='export', id=prop.type_id)[0].one().name: prop.value for prop in
                        self.db.query(Featureprop).filter_by(feature_id=feature.feature_id).all()
                    }
                    record_features.append(
                        SeqFeature(
                            FeatureLocation(featureloc.fmin, featureloc.fmax),
                            id=feature.uniquename,
                            type=CvtermService.get_query(purpose='export', id=feature.type_id),
                            strand=featureloc.strand,
                            qualifiers=qualifiers
                        )
                    )

            record = SeqRecord(seq,
                               id=org.common_name or org.infraspecific_name or org.species,
                               name=org.common_name or org.infraspecific_name or org.species,
                               description="%s %s" % (org.genus, org.species),
                               annotations={"molecule_type": "DNA"},)
            record.features = record_features
            seqs.append(record)

        return seqs

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
