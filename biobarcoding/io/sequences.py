import os.path
from typing import Tuple

from . import batch_iterator
from ..db_models import DBSessionChado, DBSession
from ..db_models.bioinformatics import Sequence, Specimen
from ..db_models.chado import Organism, Stock, Feature, StockFeature
from ..db_models.metadata import Taxon
from ..services import get_bioformat, log_exception
from ..services.bio.meta.ontologies import get_type_id
from ..services.bio.bos.sequences import Service
seq_service = Service()

BATCH_SIZE = 500
ORG_BATCH = {}
IND_BATCH = {}
ind_type_id = get_type_id(type='stock')
seq_type_id = get_type_id(type='sequence')


def split_org_name(name: str) -> Tuple[str, str, str]:
    # split a full name
    genus = species = ssp = ''
    n_split = name.split()
    if len(n_split) > 0:
        species = n_split[1-len(n_split)]
        genus = n_split[0] if n_split[0] != species else ''
        ssp = ' '.join(n_split[2:])
    return genus.strip(), species.strip(), ssp.strip()


def get_or_create(session, model, **params):
    instance = session.query(model).filter_by(**params).first()
    if not instance:
        instance = model(**params)
    return instance


# TODO: use bulk_save_objects ?
def import_file(infile, format=None, **kwargs):
    content, count = [], 0
    format = get_bioformat(infile, format)
    try:
        from ..services import seqs_parser
        stock_binds, db_taxa, db_specs, db_seqs = [], [], [], []
        for i, batch in enumerate(batch_iterator(seqs_parser(infile, format), BATCH_SIZE)):

            # to db as independent batch
            chado_seqs, chado_orgs, chado_inds = [], [], []
            for s in batch:
                org = s.annotations.get('organism')
                if org not in ORG_BATCH:
                    org_params = {}
                    org_params['genus'], org_params['species'], org_params['infraspecific_name'] = split_org_name(org)
                    ORG_BATCH[org] = get_or_create(DBSessionChado, Organism, **org_params)
                    chado_orgs.append(ORG_BATCH[org])
                    db_taxa.append(get_or_create(DBSession, Taxon,
                                                 # native_id=stock.stock_id,
                                                 # native_table='stock',
                                                 name=org))

                ind = s.name
                if ind not in IND_BATCH:
                    IND_BATCH[ind] = get_or_create(DBSessionChado, Stock, uniquename=ind, type_id=ind_type_id)
                    chado_inds.append(IND_BATCH[ind])
                    db_specs.append(get_or_create(DBSession, Specimen,
                                                  # native_id=stock.stock_id,
                                                  # native_table='stock',
                                                  name=ind))

            print('ORGS:', len(chado_orgs))
            DBSessionChado.add_all(chado_orgs)
            print('INDVS:', len(chado_inds))
            DBSessionChado.add_all(chado_inds)
            DBSessionChado.flush()

            # to db the dependent of organism part
            for s in batch:
                org = s.annotations.get('organism')
                chado_seqs.append(Feature(uniquename=s.id, name=s.name or s.description,
                                          organism_id=ORG_BATCH.get(org).organism_id,
                                          residues=str(s.seq), seqlen=len(s.seq),
                                          type_id=seq_type_id))

            DBSessionChado.add_all(chado_seqs)
            DBSessionChado.flush()

            # to db the dependent of feature part
            for s in chado_seqs:
                stock_binds.append(get_or_create(DBSessionChado, StockFeature,
                                                 stock_id=IND_BATCH.get(s.name).stock_id,
                                                 feature_id=s.feature_id,
                                                 type_id=s.type_id))

                db_seqs.append(get_or_create(DBSession, Sequence,
                                              # specimen_id=fos_specimen.id,
                                              native_id=s.feature_id,
                                              native_table='feature',
                                              name=s.uniquename))
            print('BATCH DONE')

        DBSessionChado.add_all(stock_binds)
        DBSession.add_all(db_taxa)
        DBSession.add_all(db_specs)
        DBSession.add_all(db_seqs)
    except Exception as e:
        log_exception(e)
        raise Exception(f'IMPORT sequences: file {os.path.basename(infile)} could not be imported.')
    return content, count

