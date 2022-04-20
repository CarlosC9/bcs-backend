from typing import Tuple, List

from . import MetaService
from ...main import get_orm
from ....db_models import DBSession, DBSessionChado
from ....db_models.metadata import Taxon
from ....rest import filter_parse


##
# MISCELLANEOUS
##

def get_taxonomic_ranks(rank):
    from .ontologies import CvtermService
    return CvtermService().read(cv='taxonomic_rank', cvterm=rank)[0]


def split_org_name(name: str) -> Tuple[str, str, str]:
    # split a full name
    genus = species = ssp = ''
    n_split = name.split()
    if len(n_split) > 0:
        species = n_split[1-len(n_split)]
        genus = n_split[0] if n_split[0] != species else ''
        ssp = ' '.join(n_split[2:])
    return genus.strip(), species.strip(), ssp.strip()


def build_org_name(**values):
    org = values.get('species') or ''
    if ' ' not in org:
        genus = values.get('genus', '')
        if genus == org:
            genus = ''
        ssp = values.get('infraspecific_name', '')
        org = ' '.join([genus, org, ssp])
    return org


def get_org_lineage(*org) -> list:
    from ...species_names import get_canonical_species_lineages
    return [t['canonicalName'] for t in get_canonical_species_lineages(DBSession, org)[0]]


##
# ORGANISM SERVICE
##

class Service(MetaService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('organisms')

    def get_org_name(self, organism_id=None, **values) -> str:
        # figure out and build the full name
        if organism_id:
            return self.db.query(self.orm.name).filter_by(organism_id=organism_id).one()[0].strip()
        org = values.get('organism') or values.get('name')
        if not org:
            org = build_org_name(**values)
        return org.strip()

    ##
    # CREATE
    ##

    def prepare_values(self, **values) -> dict:

        genus, species, ssp = split_org_name(
            values.get('organism') or values.get('name') or values.get('species') or '')
        values['genus'] = values.get('genus') or genus
        values['species'] = values.get('species') or species
        values['infraspecific_name'] = values.get('infraspecific_name') \
                                       or values.get('ssp./var.') or values.get('ssp') or ssp

        values['type_id'] = values.get('type_id') or values.get('cvterm_id')
        _type = values.get('type') or values.get('rank') or values.get('cvterm')
        if _type and not values.get('type_id'):
            values['type_id'] = get_taxonomic_ranks(_type)[0].cvterm_id

        return super(Service, self).prepare_values(**values)

    def check_values(self, **values) -> dict:
        values = self.check_with_gbif(**values)

        if not values.get('genus'):
            values['genus'] = ''
            # raise Exception('Missing the field "genus" for organisms.')
        if not values.get('species'):
            raise Exception('Missing the field "species" for organisms.')

        return super(Service, self).check_values(**values)

    def check_with_gbif(self, **values):
        org = self.get_org_name(**values)

        from ...species_names import get_canonical_species_info
        gbif = get_canonical_species_info(DBSession, [org])[0]
        if gbif:
            rank = gbif.get('rank', 'unknown').lower()
            values['species'] = values.get('species') or gbif.get(rank)
            # TODO: if rank is genus or above ? replace genus by rank or leave it
            # if not values.get('genus') or rank == 'genus':
            #     values['genus'] = rank
            if not values.get('infraspecific_name') and gbif.get('scientificName'):
                values['infraspecific_name'] = ' '.join(gbif.get('scientificName').split()[2:]).strip()
            try:
                values['type_id'] = values.get('type_id') or get_taxonomic_ranks(rank)[0].cvterm_id
            except:
                pass

            # gbif.update(values)
            # values = gbif

        return values

    def after_create(self, new_object, **kwargs):

        # org to bcs
        from ... import get_or_create
        try:
            app_org = get_or_create(DBSession, Taxon, name=self.get_org_name(organism_id=new_object.organism_id))
        except:
            pass

        from .taxonomies import insert_taxon
        insert_taxon(organism_id=new_object.organism_id, **kwargs)

        return new_object

    ##
    # READ
    ##

    def attach_data(self, content):
        new = super(Service, self).attach_data(content)

        if new:
            try:
                name = " ".join([new['genus'], new['species'], new['infraspecific_name'] or '']).strip()
                from ... import force_underscored
                from ...species_names import get_canonical_species_names
                new['name'] = name
                new['canonical_name'] = get_canonical_species_names(DBSession, [name])[0] or name
                new['canonical_underscored_name'] = get_canonical_species_names(DBSession, [name], underscores=True)[0] \
                                                    or force_underscored(name)
                # new['canonical_lineage'] = get_org_lineage(name)
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

        name = filter.get('organism') or filter.get('name')
        if name:
            genus, species, ssp = split_org_name(name)
            clauses += [self.orm.genus == genus, self.orm.species == species, self.orm.infraspecific_name == ssp]

        # TODO: ask by lineage: ancestor/predecessor or descendant

        if filter.get('feature_id'):
            from ....db_models.chado import Feature
            _organism_ids = self.db.query(Feature.organism_id)\
                .filter(filter_parse(Feature, [{'feature_id': filter.get('feature_id')}])).subquery()
            clauses.append(self.orm.organism_id.in_(_organism_ids))

        if filter.get('analysis_id'):
            from ....db_models.chado import AnalysisFeature, Feature
            _feature_ids = self.db.query(AnalysisFeature.feature_id)\
                .filter(filter_parse(AnalysisFeature, [{'analysis_id': filter.get('analysis_id')}])).subquery()
            _organism_ids = self.db.query(Feature.organism_id)\
                .filter(Feature.feature_id.in_(_feature_ids)).subquery()
            clauses.append(self.orm.organism_id.in_(_organism_ids))

        if filter.get('phylonode_id'):
            from ....db_models.chado import PhylonodeOrganism
            _organism_ids = self.db.query(PhylonodeOrganism.organism_id)\
                .filter(PhylonodeOrganism.phylonode_id == filter.get('phylonode_id')).subquery()
            clauses.append(self.orm.organism_id.in_(_organism_ids))

        if filter.get('phylotree_id'):
            from ....db_models.chado import Phylonode, PhylonodeOrganism
            _phylonode_ids = self.db.query(Phylonode.phylonode_id)\
                .filter(filter_parse(Phylonode, [{'phylotree_id': filter.get('phylotree_id')}])).subquery()
            _organism_ids = self.db.query(PhylonodeOrganism.organism_id)\
                .filter(PhylonodeOrganism.phylonode_id.in_(_phylonode_ids)).subquery()
            clauses.append(self.orm.organism_id.in_(_organism_ids))

        if filter.get('rank'):
            # TODO: organism.type_id ?
            from ....db_models.chado import Phylonode, PhylonodeOrganism, Cv, Cvterm
            _rank_cvterm_ids = self.db.query(Cvterm.cvterm_id)\
                .join(Cv).filter(Cv.name=='taxonomy')\
                .filter(filter_parse(Cvterm, [{'name':filter.get('rank')}])).subquery()
            _phylonode_ids = self.db.query(Phylonode.phylonode_id)\
                .filter(Phylonode.type_id.in_(_rank_cvterm_ids)).subquery()
            _ids = self.db.query(PhylonodeOrganism.organism_id)\
                .filter(PhylonodeOrganism.phylonode_id.in_(_phylonode_ids)).subquery()
            clauses.append(self.orm.organism_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(filter)
