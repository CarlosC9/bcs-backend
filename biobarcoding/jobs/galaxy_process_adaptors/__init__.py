import abc
from abc import ABC

from biobarcoding.jobs.process_adaptor import ProcessAdaptor


class GalaxyProcessAdaptor(ProcessAdaptor, ABC):
    '''The methods of the subinterfaces are all private'''

    @abc.abstractmethod
    def _complete_inputs_with_labels(self,job_context):
        raise NotImplementedError
    @abc.abstractmethod
    def _complete_with_outputs_files(self,job_context):
        raise NotImplementedError

    @abc.abstractmethod
    def _complete_with_worlflow_name(self):
        raise NotImplementedError

    def adapt_job_context(self, job_context):
        job_context= self._complete_inputs_with_labels(job_context)
        job_context['results'] = self._complete_with_outputs_files(job_context)
        job_context['process']['workflow_name'] = self._complete_with_worlflow_name()
        return job_context