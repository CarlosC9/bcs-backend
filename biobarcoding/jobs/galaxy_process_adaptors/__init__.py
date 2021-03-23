import abc
from abc import ABC

from biobarcoding.jobs.process_adaptor import ProcessAdaptor


class GalaxyProcessAdaptor(ProcessAdaptor, ABC):
    '''The methods of the subinterfaces are all private'''

    @abc.abstractmethod
    def __complete_inputs_with_labels(self,job_context):
        raise NotImplementedError
    @abc.abstractmethod
    def __complete_with_outputs_files(self,job_context):
        raise NotImplementedError

    def __complete_with_errors_files(self, job_context):
        raise NotImplementedError

    def adapt_job_context(self, job_context):
        job_context['process']['inputs']['data'] = self.__complete_inputs_with_labels(job_context)
        job_context['results'] = self.__complete_with_outputs_files(job_context)
        return job_context