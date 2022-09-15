import os.path

from Bio.SeqRecord import SeqRecord
from sqlalchemy import case, func

from . import batch_iterator
from ..db_models import DBSessionChado, DBSession
from ..db_models.bioinformatics import Sequence, Specimen
from ..db_models.chado import Organism, Stock, Feature, StockFeature, Featureloc, AnalysisFeature
from ..db_models.metadata import Taxon
from ..services import log_exception, get_encoding, tsv_pd_parser, csv_pd_parser, get_filtering, get_or_create, listify
from ..services.bio.meta.ontologies import get_type_id
from ..services.bio.meta.organisms import split_org_name, get_taxonomic_ranks

BATCH_SIZE = 500
ORG_ENTRIES, STOCK_ENTRIES, SPECIMEN_ENTRIES, SEQ_ENTRIES, ANN_ENTRIES = {}, {}, {}, {}, {}


def clear_entries():
    for _ in (ORG_ENTRIES, STOCK_ENTRIES, SPECIMEN_ENTRIES, SEQ_ENTRIES, ANN_ENTRIES):
        _.clear()


##
# SEQUENCE FILE PARSERS
##

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


##
# ANNOTATION TOOLS
##

def set_annotation(seq):

    ann = []

    # bibtex
    for ref in seq.annotations.pop('references', []):
        # TODO figure out bibtex type ?
        form_values = ref.__dict__
        loc = form_values.pop('location', None)     # TODO what for location ?
        form_values['author'] = form_values.get('author', form_values.pop('authors', None))
        ann.append({'template': {'name': 'article', 'standard': 'BibTex'},
                    'value': dict([(k, v) for k, v in form_values.items() if v])})

    # annotations
    if seq.annotations:
        ann.append({'template': {'name': 'annotations', 'standard': 'features'},
                    'value': dict([(k, v) for k, v in seq.annotations.items() if v])})

    # features
    if seq.features:
        for feat in seq.features:
            ann.append({'template': {'name': 'feature', 'standard': 'features'}, 'value': feat})

    return ann


def add_ann_entry(ann, *seq_ids):
    from .annotations import ann_value_dump
    _hash = ann_value_dump(ann)
    if ANN_ENTRIES.get(_hash):
        ANN_ENTRIES[_hash].union(seq_ids)
    else:
        ANN_ENTRIES[_hash] = {*seq_ids}


##
# SEQUENCE TOOLS
##

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


##
# MAIN PROGRAM
##

def import_file(infile, _format=None, data=None, analysis_id=None, update=False, only_residues=False, **kwargs):
    # TODO:
    #  use kwargs (p.e. organism_id, annotation_id, source)
    #  batch queries ?
    #  use bulk_save_objects ?
    try:

        clear_entries()
        global ORG_ENTRIES, STOCK_ENTRIES, SPECIMEN_ENTRIES, SEQ_ENTRIES, ANN_ENTRIES
        taxa_rows, stock_rl_rows, seq_rows, ansis_rl_rows, ansis_src_rows = [], [], [], [], []

        ind_type_id = get_type_id(type='stock')
        seq_type_id = get_type_id(type='sequence', subtype='aligned' if analysis_id else None)
        org_no_rank_id = get_taxonomic_ranks('no_rank')
        unknown_org = get_or_create(DBSessionChado, Organism, no_flush=True, genus='unknown', species='organism')
        if not unknown_org.type_id:
            unknown_org.type_id = org_no_rank_id
        DBSessionChado.add(unknown_org)

        gene = get_filtering('gene', kwargs)
        if gene:
            print(f'Extracting %s regions')
            from genbank_sequences import GenbankSeqsTools
            infile = GenbankSeqsTools.split_seqs(infile, infile + '.sp', listify(gene), _format)

        def get_seq_uniquename(s: SeqRecord):
            return f'{s.id}.a{analysis_id}' if analysis_id else s.id     # TODO: add gene ?

        try:
            _seqs = iter(data) if data else None
        except TypeError as e:
            _seqs = data if data and _format != 'fasta' else None
        for i, batch in enumerate(batch_iterator(_seqs or seqs_parser(infile, _format), BATCH_SIZE)):

            print(f'>> NEW BATCH {i+1}: {len(batch)}')
            feat_src_batch = {}
            org_batch, stock_batch, specimen_batch, feat_batch = [], [], [], []
            seq_uniquenames = [get_seq_uniquename(s) for s in batch]

            if update and only_residues:
                _new_residues, _new_seqlen = {}, {}
                for s in batch:
                    _id = get_seq_uniquename(s)
                    _new_residues[_id] = str(s.seq)
                    _new_seqlen[_id] = len(s.seq)
                    if analysis_id:
                        feat_src_batch[_id] = s.id
                DBSessionChado.query(Feature).filter(Feature.uniquename.in_(seq_uniquenames),
                                                     Feature.type_id == seq_type_id) \
                    .update({
                        Feature.residues: case(
                            _new_residues,
                            value=Feature.uniquename
                        ),
                        Feature.seqlen: case(
                            _new_seqlen,
                            value=Feature.uniquename
                        ),
                        Feature.timelastmodified: func.now()
                    }, synchronize_session=False)
                continue

            # # BATCH_READING 1: the rows that are independent or required by others
            for seq in batch:
                # STEP 1: taxa
                _org = seq.annotations.get('organism')
                if _org and _org in ORG_ENTRIES:
                    continue
                elif _org:
                    taxa_rows.append(get_or_create(DBSession, Taxon, no_flush=True, name=_org))
                    _ = dict(zip(('genus', 'species', 'infraspecific_name'), split_org_name(_org)))
                    _row = get_or_create(DBSessionChado, Organism, no_flush=True, **_)
                    if not _row.type_id:
                        _row.type_id = org_no_rank_id
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
                        _row = unknown_org
                ORG_ENTRIES[_org] = _row

            # # BATCH_READING 2: the rows that require previous data
            # STEP 0: add the required rows to the session and flush
            DBSessionChado.add_all(org_batch)
            print('ORGANISMS:', len(org_batch))
            DBSessionChado.flush()
            print('flush chado')

            # STEP 0.1: ask for and update the existent rows if updating
            if update:
                _new_name, _new_organism_id, _new_type_id, _new_is_analysis, _new_residues, _new_seqlen = \
                    {}, {}, {}, {}, {}, {}
                for s in batch:
                    _id = get_seq_uniquename(s)
                    _new_name[_id] = s.name if s.name and s.name != '<unknown name>' \
                        else s.description or s.id
                    _new_organism_id[_id] = ORG_ENTRIES[s.annotations.get('organism')].organism_id
                    _new_type_id[_id] = seq_type_id
                    _new_is_analysis[_id] = bool(analysis_id)
                    _new_residues[_id] = str(s.seq)
                    _new_seqlen[_id] = len(s.seq)
                    if analysis_id:
                        feat_src_batch[_id] = s.id
                _q = DBSessionChado.query(Feature).filter(Feature.uniquename.in_(seq_uniquenames),
                                                          Feature.type_id == seq_type_id)
                feat_batch.extend(_q.all())
                _q.update({
                    Feature.name: case(
                        _new_name,
                        value=Feature.uniquename
                    ),
                    Feature.organism_id: case(
                        _new_organism_id,
                        value=Feature.uniquename
                    ),
                    Feature.type_id: case(
                        _new_type_id,
                        value=Feature.uniquename
                    ),
                    Feature.is_analysis: case(
                        _new_is_analysis,
                        value=Feature.uniquename
                    ),
                    Feature.residues: case(
                        _new_residues,
                        value=Feature.uniquename
                    ),
                    Feature.seqlen: case(
                        _new_seqlen,
                        value=Feature.uniquename
                    ),
                    Feature.timelastmodified: func.now()
                }, synchronize_session=False)

            # STEP 0.2: ask for the existent stocks (chado individuals)
            _q_ind = [seq_name2ind(s.id)[0] for s in batch]
            _q_stock = DBSessionChado.query(Stock).filter(
                Stock.uniquename.in_([_ for _ in _q_ind if _ not in STOCK_ENTRIES]),
                Stock.type_id == ind_type_id)
            for _s in _q_stock.all():
                STOCK_ENTRIES[_s.uniquename] = _s

            _q_stock = DBSession.query(Specimen).filter(
                Specimen.name.in_([_ for _ in _q_ind if _ not in STOCK_ENTRIES]))
            for _s in _q_stock.all():
                SPECIMEN_ENTRIES[_s.name] = _s

            for seq in batch:
                # STEP 1: stock/specimen
                _org = seq.annotations.get('organism')
                _org_id = ORG_ENTRIES[_org].organism_id
                _ind, _ind_v = seq_name2ind(seq.id)
                if _ind not in STOCK_ENTRIES:
                    _ = get_or_create(DBSessionChado, Stock, no_flush=True, uniquename=_ind, type_id=ind_type_id)
                    _.organism_id = _org_id
                    STOCK_ENTRIES[_ind] = _
                    stock_batch.append(_)
                if _ind not in SPECIMEN_ENTRIES:
                    SPECIMEN_ENTRIES[_ind] = get_or_create(DBSession, Specimen, no_flush=True, name=_ind)
                    specimen_batch.append(SPECIMEN_ENTRIES[_ind])
                # STEP 2: feature (sequence)
                _id = get_seq_uniquename(seq)
                if analysis_id:
                    feat_src_batch[_id] = seq.id
                if not update:
                    _name = seq.name if seq.name and seq.name != '<unknown name>' \
                        else seq.description or seq.id
                    feat_batch.append(Feature(uniquename=_id, name=_name, organism_id=_org_id,
                                              type_id=seq_type_id, is_analysis=bool(analysis_id),
                                              residues=str(seq.seq), seqlen=len(seq.seq)))
                # STEP 3: annotations (gene/region)
                for f in seq.features:
                    if f.type == 'gene':
                        # TODO: split by gene range ?
                        for ann in f.qualifiers.get('gene'):
                            _ = {'field': {'name': 'gene', 'standard': 'NEXTGENDEM'},
                                 'value': ann.lower() if isinstance(ann, str) else ann}
                            add_ann_entry(_, seq.id)
                for ann in set_annotation(seq):
                    add_ann_entry(ann, seq.id)

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

            # STEP 0.1: ask for the existent rows if updating
            if update:
                _q = DBSession.query(Sequence.name, Sequence).filter(
                    Sequence.native_id.in_([_.feature_id for _ in feat_batch])).all()
                SEQ_ENTRIES = dict(_q + list(SEQ_ENTRIES.items()))
                seq_rows.extend(dict(_q).values())

            # STEP 0.2: ask for the source if linking any analysis
            if analysis_id:
                _srcs = DBSessionChado.query(Feature).filter(Feature.uniquename.in_(
                    [feat_src_batch.get(_f.uniquename) for _f in feat_batch])).all()
                _anl_seq_src = dict(zip([_f.uniquename for _f in feat_batch], _srcs))

            for feature in feat_batch:
                _id = feature.uniquename
                _ind, _ind_v = seq_name2ind(_id)
                # STEP 1: feature src and analysis rl
                if analysis_id:
                    try:
                        ansis_src_rows.append(get_or_create(DBSessionChado, Featureloc, no_flush=True,
                                                            feature_id=feature.feature_id,
                                                            srcfeature_id=_anl_seq_src.get(_id).feature_id))
                    except Exception as e:
                        print('Warning: Feature source could not be found.')
                    ansis_rl_rows.append(get_or_create(DBSessionChado, AnalysisFeature, no_flush=True,
                                                       analysis_id=analysis_id, feature_id=feature.feature_id))
                # STEP 2: sequence (sysadmin)
                _stock_id = STOCK_ENTRIES[_ind].stock_id
                _spec = SPECIMEN_ENTRIES[_ind]
                _spec.native_id = _stock_id
                if not update:
                    SEQ_ENTRIES[_id] = Sequence(name=_id, specimen_id=_spec.id, native_id=feature.feature_id)
                    seq_rows.append(SEQ_ENTRIES[_id])
                if not SEQ_ENTRIES.get(_id):
                    SEQ_ENTRIES[_id] = get_or_create(DBSession, Sequence,
                                                     name=_id, specimen_id=_spec.id, native_id=feature.feature_id)
                    seq_rows.append(SEQ_ENTRIES[_id])
                # STEP 3: stock relationship (individual)
                stock_rl_rows.append(get_or_create(DBSessionChado, StockFeature, no_flush=True, stock_id=_stock_id,
                                                   feature_id=feature.feature_id, type_id=feature.type_id))

        print('BATCHES DONE')
        DBSession.add_all(seq_rows)
        print('SEQUENCES:', len(seq_rows))
        DBSession.flush()
        print('flush ngd')
        DBSession.add_all(taxa_rows)
        print(f'TAXA: {len(taxa_rows)}')
        DBSessionChado.add_all(stock_rl_rows)
        print(f'STOCK_RL: {len(stock_rl_rows)}')
        DBSessionChado.add_all(ansis_rl_rows)
        print(f'ANLYS_RL: {len(ansis_rl_rows)}')
        DBSessionChado.add_all(ansis_src_rows)
        print(f'ANLYS_SRC: {len(ansis_src_rows)}')

        # TODO: run a celery task for it ?
        try:
            for k, v in ANN_ENTRIES.items():
                _list = []
                for uniquename in list(v):
                    _list.append(str(SEQ_ENTRIES[uniquename].uuid))
                ANN_ENTRIES[k] = _list
            source = get_filtering('source', kwargs)
            if source:
                _ = {'field': {'name': 'source', 'standard': 'NEXTGENDEM'},
                     'value': source}
                add_ann_entry(_, *[str(_seq.uuid) for _seqs in SEQ_ENTRIES.values() for _seq in _seqs])
            from ..tasks.system import sa_seq_ann_task
            sa_seq_ann_task.delay(ANN_ENTRIES)
        except Exception as e:
            log_exception(e)
            print('WARNING: The sequences could not be annotated.')

        # clear_entries()
    except Exception as e:
        clear_entries()
        log_exception(e)
        DBSession.rollback()
        DBSessionChado.rollback()
        raise Exception(f'IMPORT sequences: file {os.path.basename(infile)} could not be imported.')
    return None, len(seq_rows)

