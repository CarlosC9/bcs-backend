import abc

from ..rest import tm_processes


class ProcessAdaptor(abc.ABC):
    """ The only public method of Process Adaptors """

    @abc.abstractmethod
    def adapt_job_context(self, job_context):
        raise NotImplementedError

    def parse_dict_env_variables(self, d: dict):
        env_variables = ""
        for key, value in d.items():
            env_variables += f"{key}={value} "
        return env_variables.rstrip()



class ProcessAdaptorFactory:

    def get(self, jm_type, process_id):
        process_adaptor = None
        if jm_type == "galaxy":
            process_adaptor = self._get_galaxy_process_adaptor(process_id)
        elif jm_type == "ssh":
            process_adaptor = self._get_ssh_process_adaptor(process_id)
        elif jm_type == "slurm":
            process_adaptor = self._get_slurm_process_adaptor(process_id)
        elif jm_type == "cipres":
            process_adaptor = self._get_cipres_process_adaptor(process_id)
        return process_adaptor

    def _get_ssh_process_adaptor(self, process_id):
        ssh_process_adaptor = None
        if tm_processes[process_id] == "MSA ClustalW":
            from biobarcoding.jobs.ssh_process_adaptors.SSHClustalProcessAdaptor import SSHClustalProcessAdaptor
            ssh_process_adaptor = SSHClustalProcessAdaptor()
        elif tm_processes[process_id] == "PAUP Parsimony":
            from biobarcoding.jobs.ssh_process_adaptors.SSHPaupParsimonyProcessAdaptor import \
                SSHPaupParsimonyProcessAdaptor
            ssh_process_adaptor = SSHPaupParsimonyProcessAdaptor()
        elif tm_processes[process_id] == "MSA ClustalW + PAUP Parsimony":
            from biobarcoding.jobs.ssh_process_adaptors.SSHClustal_PaupParsimonyProcessAdaptor import \
                SSHClustal_PaupParsimonyProcessAdaptor
            ssh_process_adaptor = SSHClustal_PaupParsimonyProcessAdaptor()
        elif tm_processes[process_id] == "Phylogenetic Diversity Analyzer":
            from biobarcoding.jobs.ssh_process_adaptors.SSHPdaProcessAdaptor import SSHPdaProcessAdaptor
            ssh_process_adaptor = SSHPdaProcessAdaptor()
        elif tm_processes[process_id] == "PAUP Parsimony + Phylogenetic Diversity Analyzer":
            from biobarcoding.jobs.ssh_process_adaptors.SSHPaupParsimony_PdaProcessAdaptor import \
                SSHPaupParsimony_PdaProcessAdaptor
            ssh_process_adaptor = SSHPaupParsimony_PdaProcessAdaptor()
        elif tm_processes[process_id] == "MSA ClustalW + PAUP Parsimony + Phylogenetic Diversity Analyzer":
            from biobarcoding.jobs.ssh_process_adaptors.SSHClustal_PaupParsimony_PdaProcessAdaptor import \
                SSHClustal_PaupParsimony_PdaProcessAdaptor
            ssh_process_adaptor = SSHClustal_PaupParsimony_PdaProcessAdaptor()
        elif tm_processes[process_id] == "geoprocess":
            from ..jobs.ssh_process_adaptors.GeoprocessProcessAdaptor import SSHGeoprocessProcessAdaptor
            ssh_process_adaptor = SSHGeoprocessProcessAdaptor()

        return ssh_process_adaptor

    def _get_galaxy_process_adaptor(self, process_id):
        galaxy_process_adaptor = None
        if tm_processes[process_id] == "MSA ClustalW":
            from biobarcoding.jobs.galaxy_process_adaptors.GalaxyClustalProcessAdaptor import GalaxyClustalAdaptator
            galaxy_process_adaptor = GalaxyClustalAdaptator()
        return galaxy_process_adaptor

    def _get_slurm_process_adaptor(self, process_id):
        slurm_process_adaptor = None
        if tm_processes[process_id] == "MAFFT":
            from biobarcoding.jobs.ssh_process_adaptors.slurm_process_adaptors.SlurmMafftProcessAdaptor import SlurmMAFFTProcessAdaptor
            slurm_process_adaptor = SlurmMAFFTProcessAdaptor()
        elif tm_processes[process_id] == "Mr Bayes":
            from biobarcoding.jobs.ssh_process_adaptors.slurm_process_adaptors.SlurmMrBayesProcessAdaptor import SlurmMrBayesProcessAdaptor
            slurm_process_adaptor = SlurmMrBayesProcessAdaptor()
        elif tm_processes[process_id] == "BLAST":
            from biobarcoding.jobs.ssh_process_adaptors.slurm_process_adaptors.SlurmBlastProcessAdaptor import SlurmBlastProcessAdaptor
            slurm_process_adaptor = SlurmBlastProcessAdaptor()

        return slurm_process_adaptor

    def _get_cipres_process_adaptor(self, process_id):
        cipres_process_adaptor = None
        if tm_processes[process_id] == "Mr Bayes":
            from biobarcoding.jobs.cipres_process_adaptors.CipresMrBayesProcessAdaptor import CipresMrBayesProcessAdaptor
            cipres_process_adaptor = CipresMrBayesProcessAdaptor()
        elif tm_processes[process_id] == "MAFFT":
            from biobarcoding.jobs.cipres_process_adaptors.CipresMafftProcessAdaptor import CipresMafftProcessAdaptor
            cipres_process_adaptor = CipresMafftProcessAdaptor()

        return cipres_process_adaptor
