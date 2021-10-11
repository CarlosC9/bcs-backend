"""

How to manage file Inputs and Outputs NOT in BCS, but in other places
 - INPUTS: upload previously from somewhere.
 - OUTPUTS: commit into system, download it, just delete it.

"""

# ALGORITHMS
import datetime
import uuid

from sqlalchemy import Integer, Column, String, Text, ForeignKey, JSON, DateTime
from sqlalchemy.orm import relationship, backref

from . import ORMBase, GUID
from .sysadmin import Identity
from sqlalchemy.dialects.postgresql import JSONB

prefix = "jobs_"


class AlgorithmType(ORMBase):
    """ BLAST, MSA, PhylogeneticTree, GeneticDiversity"""
    __tablename__ = f"{prefix}algorithm_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))


class Algorithm(ORMBase):
    """ Description of a specific algorithm
        A process (the subject of Job execution) can call one or more algorithms

    """
    __tablename__ = f"{prefix}algorithms"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    description = Column(Text)
    algorithm_type_id = Column(Integer, ForeignKey(AlgorithmType.id), nullable=False, primary_key=False)
    algorithm_type = relationship(AlgorithmType)


# Not used
class AlgorithmConfiguration(ORMBase):
    """ Default parameters for algorithms """
    # TODO Owner of configuration: user, group or algorithm

    __tablename__ = f"{prefix}algorithm_configurations"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    description = Column(Text)
    params = Column(JSON)
    algorithm_id = Column(Integer, ForeignKey(Algorithm.id), nullable=False, primary_key=False)
    algorithm = relationship(Algorithm)


# Not used
class ProcessBrokerTemplate(ORMBase):
    """
    A template used to ensure the process is carried out.
    Steps of the template adapt to the process being executed.
    """
    __tablename__ = f"{prefix}process_broker_templates"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
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
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    description = Column(Text)
    schema_inputs = Column(JSON)  # How to specify inputs for the Process
    schema_outputs = Column(JSON)  # How to read Process outputs
    execution = Column(JSON)  # How to call the Process in a resource independent ("generic") way


class ProcessConfiguration(ORMBase):
    """ Default parameters for processes """

    __tablename__ = f"{prefix}process_configurations"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True)
    name = Column(String(80))
    description = Column(Text)
    process_id = Column(Integer, ForeignKey(Process.id), nullable=False, primary_key=False)
    process = relationship(Process)
    params = Column(JSON)


class ProcessAlgorithm(ORMBase):
    """ To register which algorithms are used in a Process """
    __tablename__ = f"{prefix}process_algorithms"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    process_id = Column(Integer, ForeignKey(Process.id), nullable=False, primary_key=False)
    process = relationship(Process, backref=backref("algorithms", cascade="all, delete-orphan"))
    algorithm_id = Column(Integer, ForeignKey(Algorithm.id), nullable=False, primary_key=False)
    algorithm = relationship(Algorithm)


# Not used
class ProcessStep(ORMBase):
    """ To enumerate the steps in a process """
    __tablename__ = f"{prefix}process_steps"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    process_id = Column(Integer, ForeignKey(Process.id), nullable=False, primary_key=False)
    process = relationship(Process, backref=backref("steps", cascade="all, delete-orphan"))


# Not used
class ProcessorArchitecture(ORMBase):
    """ Processor and architecture """
    __tablename__ = f"{prefix}processor_architectures"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))


# Not used
class OperatingSystem(ORMBase):
    """ Operating system code """
    __tablename__ = f"{prefix}operating_systems"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))


# Not used
class ComputingType(ORMBase):
    """ Computing types: sequential, MPI, OpenMP, GPU, ... """
    __tablename__ = f"{prefix}computing_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))


class JobManagementType(ORMBase):
    """ Manager types: Galaxy, ssh, EBI web service, ... """
    __tablename__ = f"{prefix}job_mgmt_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))


class ProcessInJobManagementType(ORMBase):
    """
    Adaptations of a process needed to submit it to a job manager type
    """
    __tablename__ = f"{prefix}wfs_job_mgmt_types"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    process_id = Column(Integer, ForeignKey(Process.id), nullable=False, primary_key=False)
    process = relationship(Process, backref=backref("resource_jm_adaptations", cascade="all, delete-orphan"))
    jm_type_id = Column(Integer, ForeignKey(JobManagementType.id), nullable=False, primary_key=False)
    jm_type = relationship(JobManagementType, backref=backref("process_adaptations", cascade="all, delete-orphan"))
    input_mapper = Column(JSON)
    output_mapper = Column(JSON)


class ComputeResource(ORMBase):
    """ Specific compute resource """

    __tablename__ = f"{prefix}compute_resources"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(80))
    proc_arch_id = Column(Integer, ForeignKey(ProcessorArchitecture.id), nullable=True, primary_key=False)
    proc_arch = relationship(ProcessorArchitecture)
    os_id = Column(Integer, ForeignKey(OperatingSystem.id), nullable=True, primary_key=False)
    operating_system = relationship(OperatingSystem)
    jm_type_id = Column(Integer, ForeignKey(JobManagementType.id), nullable=False, primary_key=False)
    jm_type = relationship(JobManagementType)
    jm_location = Column(JSON, nullable=False)
    jm_credentials = Column(JSON, nullable=False)
    # E.g.: jm_params = {"max_running_jobs": 4}
    #   if set, a Job is contained until the resource has number of running Jobs lower than this param
    # ALTER TABLE jobs_compute_resources ADD COLUMN jm_params JSONB;
    jm_params = Column(JSONB, nullable=True)
    agreement = Column(JSON)
    consumption_counters = Column(JSON)


class ProcessInComputeResource(ORMBase):
    """
    Specific adaptation of a Process to a ComputeResource
    """

    __tablename__ = f"{prefix}wfs_compute_resources"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    process_id = Column(Integer, ForeignKey(Process.id), nullable=False, primary_key=False)
    process = relationship(Process, backref=backref("resource_adaptations", cascade="all, delete-orphan"))
    resource_id = Column(Integer, ForeignKey(ComputeResource.id), nullable=False, primary_key=False)
    resource = relationship(ComputeResource, backref=backref("process_adaptations", cascade="all, delete-orphan"))
    native_process_id = Column(String(80))


class JobStatus(ORMBase):
    """ created, preparing_workspace, transferring_data_to_resource, submitted, completed_successfully,
        transferring_data_from_resource, cleaning_up_workspace, cancelled, paused, completed_error
    """
    __tablename__ = f"{prefix}job_statuses"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
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
    uuid = Column(GUID, unique=True, default=uuid.uuid4)

    process_id = Column(Integer, ForeignKey(Process.id), nullable=False, primary_key=False)
    process = relationship(Process)
    resource_id = Column(Integer, ForeignKey(ComputeResource.id), nullable=True, primary_key=False)
    resource = relationship(ComputeResource)
    identity_id = Column(Integer, ForeignKey(Identity.id), nullable=False, primary_key=False)
    identity = relationship(Identity, backref=backref("jobs", cascade="all, delete-orphan"))
    inputs = Column(JSON)
    status = Column(String(80))
    log = Column(Text)
    outputs = Column(JSON)


# Not used
class JobStatusLog(ORMBase):
    __tablename__ = f"{prefix}jobs_status_log"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)

    job_id = Column(Integer, ForeignKey(Job.id), nullable=False, primary_key=False)
    job = relationship(Job)
    logged_time = Column(DateTime, default=datetime.datetime.utcnow())
    logged_status_id = Column(Integer, ForeignKey(JobStatus.id), nullable=False, primary_key=False)
    logged_status = relationship(JobStatus, backref=backref("statuses", cascade="all, delete-orphan"))
