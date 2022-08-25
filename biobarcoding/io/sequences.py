import json
import os.path

from Bio.SeqRecord import SeqRecord

from . import batch_iterator
from ..db_models import DBSessionChado, DBSession, ObjectType
from ..db_models.bioinformatics import Sequence, Specimen
from ..db_models.chado import Organism, Stock, Feature, StockFeature, Featureloc, AnalysisFeature
from ..db_models.core import data_object_type_id
from ..db_models.metadata import Taxon
from ..db_models.sa_annotations import AnnotationFormField, AnnotationFormTemplate, AnnotationFormItemObjectType, \
    AnnotationTemplate, AnnotationField, AnnotationItemFunctionalObject, AnnotationFormTemplateField
from ..services import log_exception, get_encoding, tsv_pd_parser, csv_pd_parser, get_filtering, get_or_create, listify
from ..services.bio.meta.ontologies import get_type_id
from ..services.bio.meta.organisms import split_org_name, get_taxonomic_ranks

BATCH_SIZE = 500
ORG_ENTRIES = {}
STOCK_ENTRIES = {}
SPECIMEN_ENTRIES = {}
ANN_ENTRIES = {}
taxa_rows, stock_rl_rows, ann_rl_rows, ansis_rl_rows, ansis_src_rows = [], [], [], [], []
org_batch, stock_batch, specimen_batch, feat_batch, seq_batch, ann_batch = [], [], [], [], [], []


def clear_entries():
    ORG_ENTRIES.clear()
    STOCK_ENTRIES.clear()
    SPECIMEN_ENTRIES.clear()
    ANN_ENTRIES.clear()
    for _ in (taxa_rows, stock_rl_rows, ann_rl_rows, ansis_rl_rows, ansis_src_rows,
              org_batch, stock_batch, specimen_batch, feat_batch, ann_batch):
        _.clear()


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


def set_annotation(seq):

    def get_or_create_ann_form(name: str, standard: str = '', fields: list = []):
        form_template = get_or_create(DBSession, AnnotationFormTemplate, no_flush=True, name=name)
        if not form_template.id:
            DBSession.add(form_template)
            form_template.standard = form_template.standard or standard
            get_or_create(DBSession, AnnotationFormItemObjectType, no_flush=True,
                          form_item=form_template, object_type_id=data_object_type_id['sequence'])
            for field in fields:
                form_field = get_or_create(DBSession, AnnotationFormField, no_flush=True, name=field)
                DBSession.add(form_field)
                form_field.standard = form_field.standard or standard
                get_or_create(DBSession, AnnotationFormItemObjectType, no_flush=True,
                              form_item=form_field, object_type_id=data_object_type_id['sequence'])
                get_or_create(DBSession, AnnotationFormTemplateField, no_flush=True,
                              form_template=form_template, form_field=form_field)
        return form_template

    def encode_template_value(value: dict, standard: str = ''):
        form_value = value.copy()
        for k, v in value.items():
            form_field = get_or_create(DBSession, AnnotationFormField, no_flush=True, standard=standard, name=k)
            form_value[form_field.id] = form_value.get(form_field.id, form_value.pop(form_field.name, None))
        return form_value

    def to_json(arg):
        try:
            return json.dumps(arg)
        except:
            if isinstance(arg, (tuple, list, set)):
                return [to_json(_) for _ in arg]
            try:
                return to_json(arg.__dict__)
            except:
                return str(arg)

    def get_or_create_ann_template(value: dict, standard: str = '', **kwargs):
        # get_or_create for value(JSONB)
        form_value = encode_template_value(value, standard)
        instance = DBSession.query(AnnotationTemplate).filter_by(value=to_json(form_value), **kwargs).first()
        if not instance:
            instance = AnnotationTemplate(value=json.loads(to_json(form_value)), **kwargs)
        if instance not in ANN_ENTRIES:
            ANN_ENTRIES[instance] = instance
            ann_batch.append(instance)
        return instance

    ann = []

    # bibtex
    for ref in seq.annotations.pop('references', []):
        # TODO figure out bibtex type ?
        form_values = ref.__dict__
        loc = form_values.pop('location', None)     # TODO what for location ?
        form_values['author'] = form_values.get('author', form_values.pop('authors', None))
        form_template = get_or_create_ann_form(name='article', standard='BibTex', fields=form_values.keys())
        ann.append(get_or_create_ann_template(value=dict([(k, v) for k, v in form_values.items() if v]),
                                              form_template=form_template, standard='BibTex'))

    # annotations
    if seq.annotations:
        form_template = get_or_create_ann_form(name='annotations', standard='features', fields=seq.annotations.keys())
        ann.append(get_or_create_ann_template(value=dict([(k, v) for k, v in seq.annotations.items() if v]),
                                              form_template=form_template, standard='features'))

    # features
    if seq.features:
        form_template = get_or_create_ann_form(name='features', standard='features',
                                               fields=seq.features[0].__dict__.keys())
        for feat in seq.features:
            _ = []
            for k, v in feat.__dict__.items():
                if v:
                    _.append((k, to_json(v)))
            ann.append(get_or_create_ann_template(value=dict(_), form_template=form_template, standard='features'))

    return ann


def get_gene_field_id():
    gene_field = get_or_create(DBSession, AnnotationFormField, no_flush=True, name='gene')
    gene_field.standard = 'NEXTGENDEM' if not gene_field.standard else gene_field.standard
    if not gene_field.id:
        DBSession.add(gene_field)
        object_type = get_or_create(DBSession, ObjectType, no_flush=True, name='sequence')
        DBSession.add(object_type)
        DBSession.flush()
        ann_obj = get_or_create(DBSession, AnnotationFormItemObjectType, no_flush=True,
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


def import_file(infile, _format=None, data=None, analysis_id=None, update=False, **kwargs):
    # TODO:
    #  use kwargs (p.e. analysis_id, gene, organism_id)
    #  batch queries ?
    #  use bulk_save_objects ?
    try:

        clear_entries()
        gene = get_filtering('gene', kwargs)        # TODO import only this region ?
        if gene:
            from genbank_sequences import GenbankSeqsTools
            infile = GenbankSeqsTools.split_seqs(infile, infile + '.sp', listify(gene))
        ind_type_id = get_type_id(type='stock')
        seq_type_id = get_type_id(type='sequence', subtype='aligned' if analysis_id else None)
        gene_field_id = get_gene_field_id()
        org_no_rank_id = get_taxonomic_ranks('no_rank')
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
                        _row = get_or_create(DBSessionChado, Organism, no_flush=True, genus='unknown', species='organism')
                        if not _row.type_id:
                            _row.type_id = org_no_rank_id
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
                    _ = get_or_create(DBSessionChado, Stock, no_flush=True, uniquename=_ind, type_id=ind_type_id)
                    _.organism_id = _org_id
                    STOCK_ENTRIES[_ind] = _
                    stock_batch.append(_)
                if _ind not in SPECIMEN_ENTRIES:
                    SPECIMEN_ENTRIES[_ind] = get_or_create(DBSession, Specimen, no_flush=True, name=_ind)
                    specimen_batch.append(SPECIMEN_ENTRIES[_ind])
                # STEP 2: feature (sequence)
                _id = f'{seq.id}.a{analysis_id}' if analysis_id else seq.id     # TODO: add gene
                if analysis_id:
                    feat_src_batch[_id] = seq.id
                _name = seq.name if seq.name and seq.name != '<unknown name>' \
                    else seq.description or seq.id
                if update:
                    _feat = get_or_create(DBSessionChado, Feature, no_flush=True, uniquename=_id)
                    _feat.name = _name
                    _feat.organism_id = _org_id
                    _feat.type_id = seq_type_id
                    _feat.is_analysis = bool(analysis_id)
                    _feat.residues = str(seq.seq)
                    _feat.seqlen = len(seq.seq)
                else:
                    _feat = Feature(uniquename=_id, name=_name, organism_id=_org_id,
                                    type_id=seq_type_id, is_analysis=bool(analysis_id),
                                    residues=str(seq.seq), seqlen=len(seq.seq))
                feat_batch.append(_feat)
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
                for ann in set_annotation(seq):     # TODO: bulk annotations ?
                    if seq4ann_batch.get(ann):
                        seq4ann_batch[ann].add(seq.id)
                    else:
                        seq4ann_batch[ann] = set([seq.id])

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
                        ansis_src_rows.append(get_or_create(DBSessionChado, Featureloc, no_flush=True,
                                                            feature_id=feature.feature_id, srcfeature_id=src.feature_id))
                    except Exception as e:
                        print('Warning: Feature source could not be found.')
                    ansis_rl_rows.append(get_or_create(DBSessionChado, AnalysisFeature, no_flush=True,
                                                       analysis_id=analysis_id, feature_id=feature.feature_id))
                # STEP 2: sequence (sysadmin)
                _stock_id = STOCK_ENTRIES[_ind].stock_id
                _spec = SPECIMEN_ENTRIES[_ind]
                _spec.native_id = _stock_id
                if update:
                    seq_entries[_id] = get_or_create(DBSession, Sequence, no_flush=True, name=_id)
                    seq_entries[_id].specimen_id = _spec.id
                    seq_entries[_id].native_id = feature.feature_id
                else:
                    seq_entries[_id] = Sequence(name=_id, specimen_id=_spec.id, native_id=feature.feature_id)
                seq_batch.append(seq_entries[_id])
                # STEP 3: stock relationship (individual)
                stock_rl_rows.append(get_or_create(DBSessionChado, StockFeature, no_flush=True, stock_id=_stock_id,
                                                   feature_id=feature.feature_id, type_id=feature.type_id))

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
                    ann_rl_rows.append(get_or_create(DBSession, AnnotationItemFunctionalObject, no_flush=True,
                                                     annotation_id=ann.id, object_uuid=seq.uuid))
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
        DBSession.rollback()
        DBSessionChado.rollback()
        raise Exception(f'IMPORT sequences: file {os.path.basename(infile)} could not be imported.')
    return None, len(seq_batch)

