from . import MetaService
from ... import log_exception
from ...main import get_orm
from ....db_models import DBSession, DBSessionChado
from ....db_models.metadata import Taxon
from ....rest import filter_parse


##
# MISCELLANEOUS
##

def get_taxonomic_ranks(rank):
    from ....db_models.chado import Cv, Cvterm
    try:
        return DBSessionChado.query(Cvterm.cvterm_id).join(Cv) \
            .filter(Cv.name == 'taxonomic_rank', Cvterm.name == rank).one()
    except Exception as e:
        pass
    try:
        return DBSessionChado.query(Cvterm.cvterm_id).join(Cv) \
            .filter(Cv.name == 'taxonomic', Cvterm.name == rank).one()
    except Exception as e:
        pass
    try:
        return DBSessionChado.query(Cvterm.cvterm_id).join(Cv) \
            .filter(Cv.name == 'taxonomic_rank', Cvterm.name == 'no_rank').one()
    except Exception as e:
        pass
    return 1


def split_org_name(name: str) -> (str, str, str):
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
    return org.strip()


def get_orgs_lineages(*orgs) -> list:
    from ...species_names import get_canonical_species_lineages
    return [[t['canonicalName'] for t in la] for la in get_canonical_species_lineages(DBSession, orgs)]


##
# ORGANISM SERVICE
##

class Service(MetaService):

    def __init__(self):
        super(Service, self).__init__()
        self.orm = get_orm('organisms')

    def get_org_name(self, organism_id=None, **values) -> str:
        # figure out and build the full name
        if organism_id:
            return self.db.query(self.orm.name).filter_by(organism_id=organism_id).one()[0].strip()
        org = values.get('organism') or values.get('name')
        return org.strip() if org else build_org_name(**values)

    ##
    # CREATE
    ##

    def prepare_values(self, **values) -> dict:

        genus = species = ssp = ''
        if values.get('split_name'):
            genus, species, ssp = split_org_name(
                values.get('organism') or values.get('name') or values.get('species') or '')
        values['genus'] = values.get('genus') or genus
        values['species'] = species or values.get('species')
        values['infraspecific_name'] = values.get('infraspecific_name') \
                                       or values.get('ssp./var.') or values.get('ssp') or ssp

        values['type_id'] = values.get('type_id') or values.get('cvterm_id')
        _type = values.get('type') or values.get('rank') or values.get('cvterm')
        if _type and not values.get('type_id'):
            values['type_id'] = get_taxonomic_ranks(_type)

        return super(Service, self).prepare_values(**values)

    def check_values(self, **values) -> dict:
        # values = self.check_with_gbif(**values)   # TODO think if it's worth what it slows down

        if not values.get('genus'):
            values['genus'] = ''
            # raise Exception('Missing the field "genus" for organisms.')
        if not values.get('species'):
            raise Exception('Missing the field "species" for organisms.')
        if not values.get('type_id'):
            values['type_id'] = get_taxonomic_ranks('no_rank')

        return super(Service, self).check_values(**values)

    def check_with_gbif(self, **values):
        org = self.get_org_name(**values)

        from ...species_names import get_canonical_species_info
        gbif = get_canonical_species_info(DBSession, [org])[0]
        if gbif:
            rank = gbif.get('rank', 'unknown').lower()
            values['species'] = values.get('species') or gbif.get(rank)
            if not values.get('infraspecific_name') and gbif.get('scientificName'):
                values['infraspecific_name'] = ' '.join(gbif.get('scientificName').split()[2:]).strip()
            try:
                values['type_id'] = values.get('type_id') or get_taxonomic_ranks(rank)
            except:
                pass

            # gbif.update(values)
            # values = gbif

        return values

    def after_create(self, new_object, **values):
        values = super(Service, self).after_create(new_object, **values)

        # org to bcs
        from ... import get_or_create
        try:
            app_org = get_or_create(DBSession, Taxon, name=self.get_org_name(organism_id=new_object.organism_id))
        except:
            pass

        from .taxonomies import insert_taxon
        insert_taxon(organism_id=new_object.organism_id, **values)

        return values

    ##
    # READ
    ##

    def attach_data(self, *content) -> list:
        new = super(Service, self).attach_data(*content)        # TODO: thread canonical query if content is too big ?

        from ... import force_underscored
        from ...species_names import get_canonical_species_names
        names = [" ".join([_['genus'], _['species'], _['infraspecific_name'] or '']).strip() for _ in new]
        c_names = cu_names = lineages = [None] * len(names)
        try:
            c_names = get_canonical_species_names(DBSession, names)
        except Exception as e:
            print('Warning: Canonical species names could not be retrieved.')
            log_exception(e)
        try:
            cu_names = get_canonical_species_names(DBSession, names, underscores=True)
        except Exception as e:
            print('Warning: Canonical underscored species names could not be retrieved.')
            log_exception(e)
        # try:
        #     # TODO: must be cached
        #     lineages = get_orgs_lineages(*c_names)
        # except Exception as e:
        #     print('Warning: Species lineages could not be retrieved.')
        #     log_exception(e)

        for i, _ in enumerate(new):
            _['name'] = names[i]
            _['canonical_name'] = c_names[i] or names[i]
            _['canonical_underscored_name'] = cu_names[i] or force_underscored(names[i])
            # _['canonical_lineage'] = lineages[i] or []

        return new

    ##
    # EXPORT
    ##

    def data2file(self, data: list, outfile, format: str, **kwargs) -> int:
        from Bio import SeqIO
        return SeqIO.write(self.chado2biopy(data), outfile, format)

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

    def aux_filter(self, _filter: dict) -> list:
        clauses = []

        name = _filter.get('organism') or _filter.get('name')
        if name:
            genus, species, ssp = split_org_name(name)
            clauses += [self.orm.genus == genus, self.orm.species == species, self.orm.infraspecific_name == ssp]

        # TODO: ask by lineage: ancestor/predecessor or descendant

        if _filter.get('feature_id'):
            from ....db_models.chado import Feature
            _ids = self.db.query(Feature.organism_id)\
                .filter(filter_parse(Feature, [{'feature_id': _filter.get('feature_id')}]))
            clauses.append(self.orm.organism_id.in_(_ids))

        if _filter.get('analysis_id'):
            from ....db_models.chado import AnalysisFeature, Feature
            _ids = self.db.query(Feature.organism_id).join(AnalysisFeature) \
                .filter(filter_parse(AnalysisFeature, [{'analysis_id': _filter.get('analysis_id')}]))
            clauses.append(self.orm.organism_id.in_(_ids))

        if _filter.get('phylonode_id'):
            from ....db_models.chado import PhylonodeOrganism
            _ids = self.db.query(PhylonodeOrganism.organism_id)\
                .filter(PhylonodeOrganism.phylonode_id == _filter.get('phylonode_id'))
            clauses.append(self.orm.organism_id.in_(_ids))

        if _filter.get('phylotree_id'):
            from ....db_models.chado import Phylonode, PhylonodeOrganism
            _ids = self.db.query(PhylonodeOrganism.organism_id).join(Phylonode) \
                .filter(filter_parse(Phylonode, [{'phylotree_id': _filter.get('phylotree_id')}]))
            clauses.append(self.orm.organism_id.in_(_ids))

        if _filter.get('rank'):
            from ....db_models.chado import Phylonode, PhylonodeOrganism, Cv, Cvterm
            _org_ids = self.db.query(PhylonodeOrganism.organism_id).join(Phylonode).join(Cvterm).join(Cv) \
                .filter(Cv.name == 'taxonomy' or Cv.name == 'taxonomic_rank') \
                .filter(filter_parse(Cvterm, [{'name': _filter.get('rank')}]))
            _type_ids = self.db.query(Cvterm.cvterm_id).join(Cv) \
                .filter(Cv.name == 'taxonomy' or Cv.name == 'taxonomic_rank') \
                .filter(filter_parse(Cvterm, [{'name': _filter.get('rank')}]))
            from sqlalchemy import or_
            clauses.append(or_(self.orm.organism_id.in_(_org_ids), self.orm.type_id.in_(_type_ids)))

        return clauses + super(Service, self).aux_filter(_filter)
