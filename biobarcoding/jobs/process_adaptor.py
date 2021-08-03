import abc
from biobarcoding.rest import tm_processes


class ProcessAdaptor(abc.ABC):
    '''The only public method of Process Adaptors'''

    @abc.abstractmethod
    def adapt_job_context(self, job_context):
        raise NotImplementedError


class ProcessAdaptorFactory:

    def get(self, jm_type, process_id):
        process_adaptor = None
        if jm_type == "galaxy":
            process_adaptor = self._get_galaxy_process_adaptor(process_id)
        elif jm_type == "ssh":
            process_adaptor = self._get_ssh_process_adaptor(process_id)

        return process_adaptor

    def _get_ssh_process_adaptor(self, process_id):
        ssh_process_adaptor = None
        if tm_processes[process_id] == "MSA ClustalW":
            from biobarcoding.jobs.ssh_process_adaptors.SSHClustalProcessAdaptor import SSHClustalProcessAdaptor
            ssh_process_adaptor = SSHClustalProcessAdaptor()
        elif tm_processes[process_id] == "PAUP Parsimony":
            from biobarcoding.jobs.ssh_process_adaptors.SSHPaupParsimonyProcessAdaptor import SSHPaupParsimonyProcessAdaptor
            ssh_process_adaptor = SSHPaupParsimonyProcessAdaptor()
        elif tm_processes[process_id] == "MSA ClustalW + PAUP Parsimony":
            from biobarcoding.jobs.ssh_process_adaptors.SSHClustal_PaupParsimonyProcessAdaptor import SSHClustal_PaupParsimonyProcessAdaptor
            ssh_process_adaptor = SSHClustal_PaupParsimonyProcessAdaptor()
        elif tm_processes[process_id] == "Phylogenetic Diversity Analyzer":
            from biobarcoding.jobs.ssh_process_adaptors.SSHPdaProcessAdaptor import SSHPdaProcessAdaptor
            ssh_process_adaptor = SSHPdaProcessAdaptor()
        elif tm_processes[process_id] == "PAUP Parsimony + Phylogenetic Diversity Analyzer":
            from biobarcoding.jobs.ssh_process_adaptors.SSHPaupParsimony_PdaProcessAdaptor import SSHPaupParsimony_PdaProcessAdaptor
            ssh_process_adaptor = SSHPaupParsimony_PdaProcessAdaptor()
        elif tm_processes[process_id] == "MSA ClustalW + PAUP Parsimony + Phylogenetic Diversity Analyzer":
            from biobarcoding.jobs.ssh_process_adaptors.SSHClustal_PaupParsimony_PdaProcessAdaptor import SSHClustal_PaupParsimony_PdaProcessAdaptor
            ssh_process_adaptor = SSHClustal_PaupParsimony_PdaProcessAdaptor()

        return ssh_process_adaptor

    def _get_galaxy_process_adaptor(self, process_id):
        galaxy_process_adaptor = None
        if tm_processes[process_id] == "MSA ClustalW":
            from biobarcoding.jobs.galaxy_process_adaptors.GalaxyClustalProcessAdaptor import GalaxyClustalAdaptator
            galaxy_process_adaptor = GalaxyClustalAdaptator()
        return galaxy_process_adaptor
