import json
import os.path

from . import batch_iterator
from ..db_models import DBSessionChado, DBSession
from ..db_models.bioinformatics import Sequence, Specimen
from ..db_models.chado import Organism, Stock, Feature, StockFeature
from ..db_models.metadata import Taxon
from ..db_models.sa_annotations import AnnotationFormField, AnnotationField, AnnotationItemFunctionalObject
from ..services import get_bioformat, seqs_parser, log_exception
from ..services.bio.meta.ontologies import get_type_id
from ..services.bio.meta.organisms import split_org_name


BATCH_SIZE = 500
ORG_ENTRIES = {}
STOCK_ENTRIES = {}
SPECIMEN_ENTRIES = {}
ANN_ENTRIES = {}


def clear_entries():
    ORG_ENTRIES.clear()
    STOCK_ENTRIES.clear()
    SPECIMEN_ENTRIES.clear()
    ANN_ENTRIES.clear()


def get_or_create(session, model, **params):
    instance = session.query(model).filter_by(**params).first()
    if not instance:
        instance = model(**params)
    return instance


# def set_annotation(s):
#     # bibtex = s.annotations.pop('references')
#     # ann = [{'field': k, 'value': v} for k, v in s.annotations.items()]
#     # feat = [{'field': f.type, 'value': f.qualifiers} for f in s.features]
#     for f in s.features:
#         if f.type == 'gene':
#             form_field = get_or_create(DBSession, AnnotationFormField,
#                                        name='gene')
#             field = get_or_create(DBSession, AnnotationField,
#                                   form_field_id=form_field.id,
#                                   value=json.dumps(f.qualifiers.get('gene')))
#             return [form_field, field]
#     return []


def import_file(infile, format=None, **kwargs):
    # TODO: use bulk_save_objects ?
    try:

        ind_type_id = get_type_id(type='stock')
        seq_type_id = get_type_id(type='sequence')
        gene_field = get_or_create(DBSession, AnnotationFormField, name='gene')
        DBSession.add(gene_field)
        DBSession.flush()
        gene_field_id = gene_field.id

        taxa_rows, stock_rl_rows, ann_rl_rows = [], [], []
        for i, batch in enumerate(batch_iterator(seqs_parser(infile, get_bioformat(infile, format)), BATCH_SIZE)):

            print(f'>> NEW BATCH {i+1}: {len(batch)}')
            seq_entries, seq4ann_batch = {}, {}
            org_batch, stock_batch, specimen_batch, feat_batch, seq_batch, ann_batch = [], [], [], [], [], []

            # # BATCH_READING 1: the rows that are independent or required by chado:feature
            for seq in batch:
                # STEP 1: taxa
                org = seq.annotations.get('organism')
                if org not in ORG_ENTRIES:
                    if org:
                        taxa_rows.append(get_or_create(DBSession, Taxon, name=org))
                        org_params = dict(zip(('genus', 'species', 'infraspecific_name'), split_org_name(org)))
                        ORG_ENTRIES[org] = get_or_create(DBSessionChado, Organism, **org_params)
                    elif not org:
                        ORG_ENTRIES[org] = get_or_create(DBSessionChado, Organism, genus='unknown', species='organism')
                    org_batch.append(ORG_ENTRIES[org])
                # STEP 2: stock/specimen
                ind = seq.name
                if ind not in STOCK_ENTRIES:
                    STOCK_ENTRIES[ind] = get_or_create(DBSessionChado, Stock, uniquename=ind, type_id=ind_type_id)
                    stock_batch.append(STOCK_ENTRIES[ind])
                if ind not in SPECIMEN_ENTRIES:
                    SPECIMEN_ENTRIES[ind] = get_or_create(DBSession, Specimen, name=ind)
                    specimen_batch.append(SPECIMEN_ENTRIES[ind])

            # # BATCH_READING 2: the rows of chado:feature
            # STEP 1: add the required rows to the session and flush
            DBSessionChado.add_all(org_batch)
            print('ORGANISMS:', len(org_batch))
            DBSessionChado.add_all(stock_batch)
            print('STOCKS:', len(stock_batch))
            DBSessionChado.flush()
            print('flush chado')
            for seq in batch:
                # STEP 2: feature (sequence)
                org = seq.annotations.get('organism')
                feat_batch.append(Feature(uniquename=seq.id, name=seq.name or seq.description, type_id=seq_type_id,
                                          organism_id=ORG_ENTRIES.get(org).organism_id,
                                          residues=str(seq.seq), seqlen=len(seq.seq)))
                # STEP 3: annotations (gene/region)
                for f in seq.features:
                    if f.type == 'gene':
                        gene = f.qualifiers.get('gene')
                        ann = gene[-1] if isinstance(gene, (tuple, list, set)) else gene
                        if seq4ann_batch.get(ann):
                            seq4ann_batch[ann].append(seq.id)
                        else:
                            seq4ann_batch[ann] = [seq.id]
                        if ann not in ANN_ENTRIES:
                            # get_or_create for value(JSONB)
                            params = {'form_field_id': gene_field_id, 'value': json.dumps(ann)}
                            instance = DBSession.query(AnnotationField).filter_by(**params).first()
                            if not instance:
                                params['value'] = ann
                                instance = AnnotationField(**params)
                            ANN_ENTRIES[ann] = instance
                            ann_batch.append(ANN_ENTRIES[ann])
                # ann = set_annotation(s)     # TODO: bulk annotations ?

            # # BATCH_READING 3: the rows that require chado:feature
            # STEP 1: add the required rows to the session and flush
            DBSessionChado.add_all(feat_batch)
            print('FEATURES:', len(feat_batch))
            DBSessionChado.flush()
            print('flush chado')
            DBSession.add_all(specimen_batch)
            print('SPECIMENS:', len(specimen_batch))
            DBSession.flush()
            print('flush ngd')
            for feature in feat_batch:
                # STEP 2: sequence (sysadmin)
                seq = feature.uniquename
                seq_entries[seq] = get_or_create(DBSession, Sequence,
                                                 specimen_id=SPECIMEN_ENTRIES.get(feature.name).id,
                                                 native_id=feature.feature_id,
                                                 native_table='feature',
                                                 name=seq)
                seq_batch.append(seq_entries[seq])
                # STEP 3: stock relationship (individual)
                stock_rl_rows.append(get_or_create(DBSessionChado, StockFeature,
                                                   stock_id=STOCK_ENTRIES.get(feature.name).stock_id,
                                                   feature_id=feature.feature_id,
                                                   type_id=feature.type_id))

            # # BATCH_READING 4: the rows that require sysadmin:sequence
            # STEP 1: add the required rows to the session and flush
            DBSession.add_all(seq_batch)
            print('SEQUENCES:', len(seq_batch))
            DBSession.add_all(ann_batch)
            print('ANNOTATIONS:', len(ann_batch))
            DBSession.flush()
            print('flush ngd')
            # for a, seqs in seq4ann_batch.items():
            #     # STEP 2: bind seqs2ann
            #     for s in seqs:
            #         ann, seq = ANN_ENTRIES.get(a), seq_entries.get(s)
            #         if ann is None or seq is None:
            #             print(f'NOPE: {s} - {a}')
            #             continue
            #         ann_rl_rows.append(get_or_create(DBSession, AnnotationItemFunctionalObject,
            #                                          annotation=ann, object=seq))
            # print('ANNOTATED:', len(seq_batch))

        print('BATCHES DONE')
        DBSession.add_all(taxa_rows)
        print(f'TAXA: {len(taxa_rows)}')
        DBSessionChado.add_all(stock_rl_rows)
        print(f'STOCK_RL: {len(stock_rl_rows)}')
        # DBSession.add_all(ann_rl_rows)
        # print(f'ANN_RL: {len(ann_rl_rows)}')
        clear_entries()
    except Exception as e:
        clear_entries()
        log_exception(e)
        raise Exception(f'IMPORT sequences: file {os.path.basename(infile)} could not be imported.')
    return None, len(seq_batch)

