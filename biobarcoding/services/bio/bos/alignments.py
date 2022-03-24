import os.path

from Bio import AlignIO

from . import BosService
from ...main import get_orm
from ....db_models import DBSession
from ....db_models import DBSessionChado
from ....db_models.chado import Organism, Feature, AnalysisFeature
from ....db_models.bioinformatics import MultipleSequenceAlignment
from ... import get_or_create, log_exception, get_bioformat
from .sequences import Service as SeqService
from ..meta.analyses import Service as AnsisService

seq_service = SeqService()
ansis_service = AnsisService()


##
# ALIGNMENT SERVICE
##
class Service(BosService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('alignments')
        self.bos = 'multiple-sequence-alignment'

    ##
    # CREATE
    ##

    def prepare_values(self, **values):
        values.update(ansis_service.check_values(**values))
        return super(Service, self).check_values(**values)


    def check_values(self, **values) -> dict:

        if not values.get('program'):
            values['program'] = 'Multiple sequence alignment'

        values.update(ansis_service.check_values(**values))
        return super(Service, self).check_values(**values)

    def after_create(self, new_object, **values):
        super(Service, self).after_create(new_object, **values)

        fos_msa = get_or_create(DBSession, MultipleSequenceAlignment,
                                native_id=new_object.analysis_id,
                                native_table='analysis',
                                name=new_object.name)

        return values

    ##
    # READ
    ##
    
    def attach_data(self, content):
        new = super(Service, self).attach_data(content)

        from sqlalchemy.sql.expression import func, distinct
        info = self.db.query(func.max(func.length(Feature.residues)),
                             func.count(AnalysisFeature.analysis_id),
                             func.array_agg(distinct(Organism.name))) \
            .select_from(AnalysisFeature).join(Feature).join(Organism) \
            .filter(AnalysisFeature.analysis_id == content['analysis_id']) \
            .group_by(AnalysisFeature.analysis_id)
        content['seqlen'], content['seqnum'], content['taxa'] = info.one()

        return new

    ##
    # DELETE
    ##

    def after_delete(self, *content, **kwargs):
        ids = [a.analysis_id for a in content]
        query = DBSession.query(MultipleSequenceAlignment).filter(MultipleSequenceAlignment.native_id.in_(ids))
        return query.count()
        # TODO: check why all rows are deleted in bcs without filtering
        # return query.delete(synchronize_session='fetch')

    ##
    # IMPORT
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

    def import_file(self, infile, format=None, **kwargs):
        content, count = None, 0
        format = get_bioformat(infile, format)
        try:
            # TODO: try every available format ('clustal', 'emboss', 'fasta', 'fasta-m10', 'ig', 'msf', 'nexus', 'phylip', 'phylip-sequential', 'phylip-relaxed', 'stockholm', 'mauve')
            # check aligned file
            content_file = AlignIO.read(infile, format or 'fasta')
            # Set missing default values
            kwargs['programversion'] = kwargs.get('programversion') or '(Imported file)'
            kwargs['sourcename'] = kwargs.get('sourcename') or os.path.basename(infile)
            # Analysis row could exist for jobs, so get or create
            try:
                unique_keys = ['job_id'] if kwargs.get('job_id') else ('program', 'programversion', 'sourcename')
                content, count = ansis_service.get_query(purpose='annotate', **{k: kwargs[k] for k in unique_keys if k in kwargs})
                content = content.one()
                # Analysis row could be created by importing the results of other jobs, and without register as msa in bcs
                self.after_create(content, **kwargs)
            except:
                content, status = self.create(**kwargs)
            # Read and import file content
            self.msafile2chado(content, content_file)
        except Exception as e:
            log_exception(e)
            raise Exception(f'IMPORT alignments: The file {os.path.basename(infile)} could not be imported.')
        return content, count

    ##
    # EXPORT
    ##

    def export_file(self, format='fasta', values={}, **kwargs):
        content, count = [], 0
        try:
            if format in ('fasta', 'nexus'):
                header = {'header': values.pop('header', None)}
                query, count = self.get_query(purpose='export', values=values, **kwargs)
                for a in query.all():
                    c, cc = seq_service.export_file(analysis_id=a.analysis_id, values=header)
                    content.append(c)
                    count += cc
            else:
                raise Exception(f'EXPORT alignments: The format {format} could not be exported.')
        except Exception as e:
            log_exception(e)
            raise Exception(f'EXPORT alignments: The alignment could not be exported.')
        return content, count
