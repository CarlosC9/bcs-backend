import os.path

from Bio import AlignIO

from . import BosService
from ...main import get_orm
from ....db_models import DBSession
from ....db_models import DBSessionChado
from ....db_models.chado import Organism, Feature, AnalysisFeature
from ....db_models.bioinformatics import MultipleSequenceAlignment
from ... import get_or_create, log_exception
from .sequences import Service as SeqService
from ..meta.analyses import Service as AnsisService

seq_service = SeqService()
# ansis_service = AnsisService()  # TODO keep an eye on the BOS and Analysis services, they might be needed


##
# ALIGNMENT SERVICE
##
class Service(BosService, AnsisService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('alignments')
        self.bos = 'multiple-sequence-alignment'
        self.formats = ['fasta', 'clustal', 'nexus', 'phylip', 'emboss', 'fasta-m10', 'ig', 'maf', 'mauve', 'msf', 'phylip-sequential', 'phylip-relaxed', 'stockholm']

    ##
    # CREATE
    ##

    def prepare_values(self, **values):
        # values.update(AnsisService.prepare_values(self, **values))
        return super(Service, self).prepare_values(**values)

    def check_values(self, **values) -> dict:

        if not values.get('program'):
            values['program'] = 'Multiple sequence alignment'

        # values.update(super(Service, self).check_values(**values))
        return super(Service, self).check_values(**values)

    def after_create(self, new_object, **values):
        values = super(Service, self).after_create(new_object, **values)

        fos_msa = get_or_create(DBSession, MultipleSequenceAlignment,
                                native_id=new_object.analysis_id,
                                native_table='analysis',
                                name=new_object.name)

        return values

    ##
    # READ
    ##
    
    def attach_data(self, *content):
        new = super(Service, self).attach_data(*content)

        _ids = [_['analysis_id'] for _ in new]
        info = [None, None, []] * len(new)
        try:
            from sqlalchemy.sql.expression import case, func, distinct
            _ord = case({_id: index for index, _id in enumerate(_ids)},
                        value=AnalysisFeature.analysis_id)  # TODO test order
            info = self.db.query(func.max(func.length(Feature.residues)),
                                 func.count(AnalysisFeature.analysis_id),
                                 func.array_agg(distinct(Organism.name))) \
                .select_from(AnalysisFeature).join(Feature).join(Organism) \
                .filter(AnalysisFeature.analysis_id.in_(_ids)) \
                .group_by(AnalysisFeature.analysis_id).order_by(_ord).all()
        except Exception as e:
            print('Error: Additional data could not be attached.')
            log_exception(e)
        for i, _ in enumerate(new):
            _['seqlen'], _['seqnum'], _['taxa'] = info[i]

        return new

    ##
    # DELETE
    ##

    def delete_related(self, *content, **kwargs):
        ids = [a.analysis_id for a in content]

        from .phylotrees import Service as PhyService
        PhyService().delete(filter={'analysis_id': ids})
        seq_service.delete(filter={'analysis_id': ids})

        query = DBSession.query(MultipleSequenceAlignment).filter(MultipleSequenceAlignment.native_id.in_(ids))
        return len([DBSession.delete(row) for row in query.all()])

    ##
    # DEPRECATED IMPORT
    ##

    def seq_org_id(self, name):
        try:
            return self.db.query(Organism).join(Feature).filter(Feature.uniquename == name).one().organism_id
        except:
            return None

    def bind2src(self, feature, srcname):
        try:
            src = seq_service.get_query(purpose='annotate', uniquename=srcname)[0].one()
            from ....db_models.chado import Featureloc
            rl = get_or_create(self.db, Featureloc, feature_id=feature.feature_id, srcfeature_id=src.feature_id)
        except:
            rl = None
        return rl

    def bind2ansis(self, msa, feature):
        try:
            rl = get_or_create(self.db, AnalysisFeature, analysis_id=msa.analysis_id, feature_id=feature.feature_id)
        except:
            rl = None
        return rl

    def msafile2chado(self, msa, seqs):
        for seq in seqs:
            feature, count = seq_service.create(
                uniquename=f'{seq.id}.msa{msa.analysis_id}',
                name=seq.id,
                stock=seq.id,
                residues=str(seq.seq),
                organism_id=self.seq_org_id(seq.id),
                type='sequence', subtype='aligned')
            self.bind2src(feature, seq.id)
            self.bind2ansis(msa, feature)
        return msa

    ##
    # IMPORT
    ##

    def read_infile(self, file, _format) -> any:
        return AlignIO.read(file, _format)    # too slow
        # from Bio import SeqIO
        # return next(SeqIO.parse(file, _format))

    def import_file(self, infile, format=None, **kwargs):
        content_file, _format = self.check_infile(infile, format)
        # Set missing default values
        kwargs['programversion'] = kwargs.get('programversion') or '(Imported file)'
        kwargs['sourcename'] = kwargs.get('sourcename') or os.path.basename(infile)
        # Analysis row could exist for jobs, so get or create
        try:
            unique_keys = ['job_id'] if kwargs.get('job_id') else ('program', 'programversion', 'sourcename')
            content, count = self.get_query(purpose='annotate', **{k: kwargs[k] for k in unique_keys if k in kwargs})
            content = content.one()
            # Analysis row could be created by importing the results of other jobs, and without register as msa in bcs
            self.after_create(content, **kwargs)
        except:
            content, count = self.create(**kwargs)
        # Read and import file content
        from ....io.sequences import import_file
        return import_file(infile, _format, data=content_file, analysis_id=content.analysis_id)

    ##
    # EXPORT
    ##

    def chado2biopy(self, aln, header: str = None) -> list:
        query, count = seq_service.get_query(purpose='export', filter={'analysis_id': aln.analysis_id})
        return seq_service.chado2biopy(query.all(), header)

    def data2file(self, alns: list, outfile, format: str, **kwargs) -> int:
        files, count = [], 0
        header = kwargs.get('values', {}).get('header')
        for a in alns:
            file = outfile
            if len(alns) > 1:
                file = f'{file}_{a.analysis_id}'
            from Bio.Align import MultipleSeqAlignment
            aln = MultipleSeqAlignment(self.chado2biopy(a, header))
            AlignIO.write(aln, file, format)
            files.append(file)
            count += 1
        if len(alns) > 1:
            from ....common.helpers import zip_files
            zip_files(outfile, files)
        return count

    ##
    # GETTER AND OTHERS
    ##

    def get_query(self, **kwargs):
        # return AnsisService.get_query(query=BosService.get_query(**kwargs)[0], **kwargs)
        return super(Service, self).get_query(**kwargs)

    # def aux_filter(self, filter):
    #     return ansis_service.aux_filter(filter) + super(Service, self).aux_filter(filter)
