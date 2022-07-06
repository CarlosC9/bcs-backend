import json
import os.path

from Bio.SeqRecord import SeqRecord

from . import batch_iterator
from ..db_models import DBSessionChado, DBSession, ObjectType
from ..db_models.bioinformatics import Sequence, Specimen
from ..db_models.chado import Organism, Stock, Feature, StockFeature, Featureloc, AnalysisFeature
from ..db_models.metadata import Taxon
from ..db_models.sa_annotations import AnnotationFormField, AnnotationField, AnnotationItemFunctionalObject, \
    AnnotationFormItemObjectType
from ..services import log_exception, get_encoding, tsv_pd_parser, csv_pd_parser
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


def gff_parser(file) -> SeqRecord:
    from BCBio import GFF
    with open(file, encoding=get_encoding(file)) as f:
        for rec in GFF.parse(f):
            yield rec


def decode_fasta_header(seq: SeqRecord) -> SeqRecord:
    pattern = '^(?P<id>(?P<name>[^.\s]+)(\.(?P<version>\S*))?)' \
              '(\s*\[organism=(?P<organism>.+)\])?(\s+(?P<features>.+);)?(\s+(?P<origin>.+))?$'
    # p.e.: >Seq1_1.2 [organism=Pinus canariensis] maturase K (matk) gene, partial sequence, partial cds; chloroplast
    import re
    meta = re.match(pattern, seq.description)
    if meta:
        meta = meta.groupdict()
        seq.name = meta.get('name') if seq.name == seq.id else seq.name
        organism = meta.get('organism')
        seq.annotations['organism'] = organism
        origin = meta.get('origin')
        seq.annotations['source'] = f'{origin} {organism}'
        features = meta.get('features')
        # TODO ? seq.molecule_type
        if features:
            if isinstance(features, str):
                features = features.split(',')
            from Bio.SeqFeature import SeqFeature
            seq.features = [SeqFeature(type=f.split()[-1],
                                       qualifiers={f.split()[-1]:f[:-len(f.split()[-1])].strip()})
                            for f in features]
    return seq


def seqs_parser(file, _format='fasta') -> SeqRecord:
    if _format == 'gff':
        return gff_parser(file)
    if _format == 'tsv':
        return SeqRecord(**tsv_pd_parser(file))
    if _format == 'csv':
        return SeqRecord(**csv_pd_parser(file))
    from Bio import SeqIO
    with open(file, encoding=get_encoding(file)) as fo:
        for rec in SeqIO.parse(fo, _format):
            if _format == 'fasta':
                rec = decode_fasta_header(rec)
            yield rec


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


def get_gene_field_id():
    gene_field = get_or_create(DBSession, AnnotationFormField, name='gene')
    DBSession.add(gene_field)
    object_type = get_or_create(DBSession, ObjectType, name='sequence')
    DBSession.add(object_type)
    DBSession.flush()
    ann_obj = get_or_create(DBSession, AnnotationFormItemObjectType,
                            form_item_id=gene_field.id,
                            object_type_id=object_type.id)
    DBSession.add(ann_obj)
    return gene_field.id


def seq_name2ind(seq: str) -> (str, str):
    _ = seq.split('.')
    return _[0], ''.join(_[1:])


def find_org(seq) -> any:  # try to find out the organism from db
    _ind, _ind_v = seq_name2ind(seq.id)
    try:
        return DBSessionChado.query(Organism).join(Feature).distinct(Organism.organism_id) \
            .filter(Feature.uniquename.like(_ind + '%')).one()
    except Exception as e:
        print('Warning: Organism information missing for %s. Stock not found.' % seq.id)
    try:
        return DBSessionChado.query(Organism).join(Stock).distinct(Organism.organism_id) \
            .filter(Stock.uniquename.like(_ind + '%')).one()
    except Exception as e:
        print('Warning: Organism information missing for %s. Related Features not found.' % seq.id)
    return None


def import_file(infile, _format=None, data=None, analysis_id=None, **kwargs):
    # TODO:
    #  use kwargs (p.e. analysis_id, gene, organism_id)
    #  batch queries ?
    #  use bulk_save_objects ?
    try:

        ind_type_id = get_type_id(type='stock')
        seq_type_id = get_type_id(type='sequence', subtype='aligned' if analysis_id else None)
        gene_field_id = get_gene_field_id()
        taxa_rows, stock_rl_rows, ann_rl_rows, ansis_rl_rows, ansis_src_rows = [], [], [], [], []

        try:
            _seqs = iter(data) if data else None
        except TypeError as e:
            _seqs = data if data and _format != 'fasta' else None
        for i, batch in enumerate(batch_iterator(_seqs or seqs_parser(infile, _format), BATCH_SIZE)):

            print(f'>> NEW BATCH {i+1}: {len(batch)}')
            seq_entries, seq4ann_batch, feat_src_batch = {}, {}, {}
            org_batch, stock_batch, specimen_batch, feat_batch, seq_batch, ann_batch = [], [], [], [], [], []

            # # BATCH_READING 1: the rows that are independent or required by others
            for seq in batch:
                # STEP 1: taxa
                _org = seq.annotations.get('organism')
                if _org and _org in ORG_ENTRIES:
                    continue
                elif _org:
                    taxa_rows.append(get_or_create(DBSession, Taxon, name=_org))
                    _ = dict(zip(('genus', 'species', 'infraspecific_name'), split_org_name(_org)))
                    _row = get_or_create(DBSessionChado, Organism, **_)
                    org_batch.append(_row)
                else:
                    _row = find_org(seq)
                    if _row and _row.organism_id:
                        seq.annotations = seq.annotations.copy()
                        seq.annotations['organism'] = _org = _row.organism_id
                    elif _org in ORG_ENTRIES:
                        print('Warning: Unknown organism for', seq.id)
                        continue
                    else:
                        print('Warning: Unknown organism for', seq.id)
                        _row = get_or_create(DBSessionChado, Organism, genus='unknown', species='organism')
                        org_batch.append(_row)
                ORG_ENTRIES[_org] = _row

            # # BATCH_READING 2: the rows that require previous data
            # STEP 0: add the required rows to the session and flush
            DBSessionChado.add_all(org_batch)
            print('ORGANISMS:', len(org_batch))
            DBSessionChado.flush()
            print('flush chado')
            for seq in batch:
                # STEP 1: stock/specimen
                _org = seq.annotations.get('organism')
                _org_id = ORG_ENTRIES[_org].organism_id
                _ind, _ind_v = seq_name2ind(seq.id)
                if _ind not in STOCK_ENTRIES:
                    _ = get_or_create(DBSessionChado, Stock, uniquename=_ind, type_id=ind_type_id)
                    _.organism_id = _org_id
                    STOCK_ENTRIES[_ind] = _
                    stock_batch.append(_)
                if _ind not in SPECIMEN_ENTRIES:
                    SPECIMEN_ENTRIES[_ind] = get_or_create(DBSession, Specimen, name=_ind)
                    specimen_batch.append(SPECIMEN_ENTRIES[_ind])
                # STEP 2: feature (sequence)
                _id = f'{seq.id}.a{analysis_id}' if analysis_id else seq.id
                if analysis_id:
                    feat_src_batch[_id] = seq.id
                _name = seq.name if seq.name and seq.name != '<unknown name>' \
                            else seq.description or seq.id
                feat_batch.append(Feature(uniquename=_id, name=_name, organism_id=_org_id,
                                          type_id=seq_type_id, is_analysis=bool(analysis_id),
                                          residues=str(seq.seq), seqlen=len(seq.seq)))
                # STEP 3: annotations (gene/region)
                for f in seq.features:
                    if f.type == 'gene':
                        # TODO: split by gene range
                        gene = f.qualifiers.get('gene')
                        ann = gene[-1] if isinstance(gene, (tuple, list, set)) else gene
                        if seq4ann_batch.get(ann):
                            seq4ann_batch[ann].add(seq.id)
                        else:
                            seq4ann_batch[ann] = set([seq.id])
                        if ann not in ANN_ENTRIES:
                            # get_or_create for value(JSONB)
                            params = {'form_field_id': gene_field_id, 'value': json.dumps(ann)}
                            instance = DBSession.query(AnnotationField).filter_by(**params).first()
                            if not instance:
                                params['value'] = ann
                                instance = AnnotationField(**params)
                            ANN_ENTRIES[ann] = instance
                            ann_batch.append(instance)
                # ann = set_annotation(s)     # TODO: bulk annotations ?

            # # BATCH_READING 3: the rows that require chado:feature
            # STEP 0: add the required rows to the session and flush
            DBSessionChado.add_all(stock_batch)
            print('STOCKS:', len(stock_batch))
            DBSessionChado.add_all(feat_batch)
            print('FEATURES:', len(feat_batch))
            DBSessionChado.flush()
            print('flush chado')
            DBSession.add_all(specimen_batch)
            print('SPECIMENS:', len(specimen_batch))
            DBSession.flush()
            print('flush ngd')
            for feature in feat_batch:
                _id = feature.uniquename
                _ind, _ind_v = seq_name2ind(_id)
                # STEP 1: feature src and analysis rl
                if analysis_id:
                    try:
                        src = DBSessionChado.query(Feature).filter(Feature.uniquename == feat_src_batch.get(_id)).one()
                        ansis_src_rows.append(Featureloc(feature_id=feature.feature_id, srcfeature_id=src.feature_id))
                    except Exception as e:
                        print('Warning: Feature source could not be found.')
                    ansis_rl_rows.append(AnalysisFeature(analysis_id=analysis_id, feature_id=feature.feature_id))
                # STEP 2: sequence (sysadmin)
                _stock_id = STOCK_ENTRIES[_ind].stock_id
                _spec = SPECIMEN_ENTRIES[_ind]
                _spec.native_id = _stock_id
                seq_entries[_id] = Sequence(name=_id, specimen_id=_spec.id, native_id=feature.feature_id)
                seq_batch.append(seq_entries[_id])
                # STEP 3: stock relationship (individual)
                stock_rl_rows.append(StockFeature(stock_id=_stock_id,
                                                  feature_id=feature.feature_id,
                                                  type_id=feature.type_id))

            # # BATCH_READING 4: the rows that require sysadmin:sequence
            # STEP 0: add the required rows to the session and flush
            DBSession.add_all(seq_batch)
            print('SEQUENCES:', len(seq_batch))
            DBSession.add_all(ann_batch)
            print('ANNOTATIONS:', len(ann_batch))
            DBSession.flush()
            print('flush ngd')
            for a, seqs in seq4ann_batch.items():
                # STEP 1: bind seqs2ann
                ann = ANN_ENTRIES.get(a)
                if not ann:
                    print(f'ANN NOT FOUND: {a}')
                    continue
                for s in seqs:
                    seq = seq_entries.get(s)
                    if not seq:
                        print(f'ANN_SEQ NOT FOUND: {s} - {a}')
                        continue
                    ann_rl_rows.append(AnnotationItemFunctionalObject(annotation_id=ann.id, object_uuid=seq.uuid))
            print('ANNOTATED:', len(seq_batch))

        print('BATCHES DONE')
        DBSession.add_all(taxa_rows)
        print(f'TAXA: {len(taxa_rows)}')
        DBSessionChado.add_all(stock_rl_rows)
        print(f'STOCK_RL: {len(stock_rl_rows)}')
        DBSession.add_all(ann_rl_rows)
        print(f'ANN_RL: {len(ann_rl_rows)}')
        DBSessionChado.add_all(ansis_rl_rows)
        print(f'ANLYS_RL: {len(ansis_rl_rows)}')
        DBSessionChado.add_all(ansis_src_rows)
        print(f'ANLYS_SRC: {len(ansis_src_rows)}')
        clear_entries()
    except Exception as e:
        clear_entries()
        log_exception(e)
        raise Exception(f'IMPORT sequences: file {os.path.basename(infile)} could not be imported.')
    return None, len(seq_batch)

