"""

How to manage file Inputs and Outputs NOT in BCS, but in other places
 - INPUTS: upload previously from somewhere.
 - OUTPUTS: commit into system, download it, just delete it.

"""

# ALGORITHMS
import datetime

from sqlalchemy import Integer, Column, String, Text, ForeignKey, JSON, DateTime
from sqlalchemy.orm import relationship, backref

from biobarcoding.db_models import ORMBase, GUID

prefix = "jobs_"


class AlgorithmType(ORMBase):
    """ BLAST, MSA, PhylogeneticTree, GeneticDiversity"""
    __tablename__ = f"{prefix}algorithm_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class Algorithm(ORMBase):
    """ Description of a specific algorithm"""
    __tablename__ = f"{prefix}algorithms"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    description = Column(Text)
    algorithm_type_id = Column(Integer, ForeignKey(AlgorithmType.id), nullable=False, primary_key=False)
    algorithm_type = relationship(AlgorithmType)


class ProcessBrokerTemplate(ORMBase):
    """
    A template used to ensure the process is carried out.
    Steps of the template adapt to the process being executed. For instance,
    """
    __tablename__ = f"{prefix}process_broker_templates"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class Process(ORMBase):
    """ Repeatable processes executed locally or in a compute resource to process bioinformatic information
      Computation
      - Multiple alignment
      - BLAST
      - Phylogenetic tree
      - Phylogenetic diversity
      Possibly long operations
      - Import/export sequences
      - Import/export multiple aligment
      - Import/export phylogenetic tree
      - Prepare complex visualization
      Maintainance
      - Database clean-up
      - Batch notifications
      - Single notification
    """
    __versioned__ = {}
    __tablename__ = f"{prefix}processes"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    description = Column(Text)
    schema_inputs = Column(JSON)  # How to specify inputs for the Process
    schema_outputs = Column(JSON)  # How to read Process outputs
    execution = Column(JSON)  # How to call the Process in a resource independent ("generic") way


class ProcessAlgorithm(ORMBase):
    """ To register which algorithms are used in a Process """
    __tablename__ = f"{prefix}process_algorithms"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    process_id = Column(Integer, ForeignKey(Process.id), nullable=False, primary_key=False)
    process = relationship(Process, backref=backref("algorithms", cascade="all, delete-orphan"))
    algorithm_id = Column(Integer, ForeignKey(Algorithm.id), nullable=False, primary_key=False)
    algorithm = relationship(Algorithm)


class ProcessStep(ORMBase):
    """ To enumerate the steps in a process """
    __tablename__ = f"{prefix}process_steps"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    process_id = Column(Integer, ForeignKey(Process.id), nullable=False, primary_key=False)
    process = relationship(Process, backref=backref("steps", cascade="all, delete-orphan"))


class ProcessorArchitecture(ORMBase):
    """ Processor and architecture """
    __tablename__ = f"{prefix}processor_architectures"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class OperatingSystem(ORMBase):
    """ Operating system code """
    __tablename__ = f"{prefix}operating_systems"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class ComputingType(ORMBase):
    """ Computing types: sequencial, MPI, OpenMP, GPU, ... """
    __tablename__ = f"{prefix}computing_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class JobManagementType(ORMBase):
    """ Manager types: Galaxy, ssh, EBI web service, ... """
    __tablename__ = f"{prefix}job_mgmt_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class ComputeResource(ORMBase):
    """ Specific compute resource """

    __tablename__ = f"{prefix}compute_resources"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    proc_arch_id = Column(Integer, ForeignKey(ProcessorArchitecture.id), nullable=True, primary_key=False)
    proc_arch = relationship(ProcessorArchitecture)
    os_id = Column(Integer, ForeignKey(OperatingSystem.id), nullable=True, primary_key=False)
    operating_system = relationship(OperatingSystem)
    jm_type_id = Column(Integer, ForeignKey(JobManagementType.id), nullable=False, primary_key=False)
    jm_type = relationship(JobManagementType)
    jm_location = Column(JSON, nullable=False)
    jm_credentials = Column(JSON, nullable=False)
    agreement = Column(JSON)
    consumption_counters = Column(JSON)


class ProcessInComputeResource(ORMBase):
    """
    Specific adaptation of a Process to a ComputeResource
    """

    __tablename__ = f"{prefix}wfs_compute_resources"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    process_id = Column(Integer, ForeignKey(Process.id), nullable=False, primary_key=False)
    process = relationship(Process, backref=backref("resource_adaptations", cascade="all, delete-orphan"))
    resource_id = Column(Integer, ForeignKey(ComputeResource.id), nullable=False, primary_key=False)
    resource = relationship(ComputeResource, backref=backref("process_adaptations", cascade="all, delete-orphan"))


class JobStatus(ORMBase):
    """ created, preparing_workspace, transferring_data_to_resource, submitted, completed_successfully,
        transferring_data_from_resource, cleaning_up_workspace, cancelled, paused, completed_error
    """
    __tablename__ = f"{prefix}job_statuses"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))


class Job(ORMBase):
    """ Core entity in the Algorithmic domain.
    Processes are launched into ComputeResources with some parameters and under an Identity
    Then, jobs can be checked, cancelled and results obtained
    Also, listing of filtered/ordered jobs may be retrieved
    NOTE: IF NO COMPUTERESOURCE IS SPECIFIED, THE JOB IS KEPT ON HOLD. LATER, A COMPUTE RESOURCE CAN BE ASSIGNED
    """
    __tablename__ = f"{prefix}jobs"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)

    process_id = Column(Integer, ForeignKey(Process.id), nullable=False, primary_key=False)
    process = relationship(Process)
    resource_id = Column(Integer, ForeignKey(ComputeResource.id), nullable=True, primary_key=False)
    resource = relationship(ComputeResource)
    status = Column(String(80))
    log = Column(Text)
    inputs = Column(JSON)
    outputs = Column(JSON)


class JobStatusLog(ORMBase):
    __tablename__ = f"{prefix}jobs_status_log"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)

    job_id = Column(Integer, ForeignKey(Job.id), nullable=False, primary_key=False)
    job = relationship(Job)
    logged_time = Column(DateTime, default=datetime.datetime.utcnow())
    logged_status_id = Column(Integer, ForeignKey(JobStatus.id), nullable=False, primary_key=False)
    logged_status = relationship(JobStatus, backref=backref("statuses", cascade="all, delete-orphan"))
