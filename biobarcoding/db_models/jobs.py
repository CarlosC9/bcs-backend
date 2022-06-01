"""

How to manage file Inputs and Outputs NOT in the App, but in other places
 - INPUTS: upload previously from somewhere.
 - OUTPUTS: commit into system, download it, just delete it.

"""

# ALGORITHMS
from enum import Enum
import uuid
import os
import json
from datetime import datetime

from sqlalchemy import Integer, Column, String, Text, ForeignKey, JSON, Boolean, DateTime, \
    UniqueConstraint, CheckConstraint
from sqlalchemy.dialects.postgresql import JSONB, ENUM
from sqlalchemy.orm import relationship, backref

from . import ORMBase, GUID
from .sysadmin import Identity
from ..common import generate_json

prefix = "jobs_"


class AlgorithmType(ORMBase):
    """ BLAST, MSA, PhylogeneticTree, GeneticDiversity, Geoprocessing, Forecasting, ... """
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


class Process(ORMBase):
    """ Repeatable processes executed locally or in a compute resource to process
        scientific / technical information
      Computation
      - Forecasting
      - Geoprocesses
      - Multiple alignment
      - BLAST
      - Phylogenetic tree
      - Phylogenetic diversity
      Possibly long operations
      - Import/export geographic layers
      - Import/export sequences
      - Import/export multiple aligment
      - Import/export phylogenetic tree
      - Prepare a complex report or visualization
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


class JobManagementType(ORMBase):
    """ Manager types: ssh, ssh-slurm, Galaxy, EBI web service, ... """
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
    jm_type_id = Column(Integer, ForeignKey(JobManagementType.id), nullable=False, primary_key=False)
    jm_type = relationship(JobManagementType)
    jm_location = Column(JSON, nullable=False)
    jm_credentials = Column(JSON, nullable=False)
    # E.g.: jm_params = {"max_running_jobs": 4}
    #   if set, a Job is contained until the resource has number of running Jobs lower than this param
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
    # JOB STATUS
    # * Initially it should be "created" (although the submission is to the initial Celery task, currently "prepare")
    # * inform users (Celery tasks change it to help tracing)
    # * If user sets the status to "cancelling", Celery tasks (definitions.py) must end themselves and
    #   return "cancelled" to jump to that task, that after a final cleanup, would change the status from
    #   "cancelling" to "cancelled".
    # * Jobs can be hard-deleted when in certain statuses: cancelled, error. If not, soft-deleted (see "deleted" below)
    #
    # TASKS ("wf1" in "definitions.py") -> JOB.STATUS (this table) -> JOB_STATUS table (see tm_job_statuses)
    # * <from POST> -> created -> created
    # * prepare -> preparing_workspace -> preparing_workspace
    # * export -> export -> exporting_to_supported_file_formats
    # * transfer_data -> transfer_data_to_repository -> transferring_data_to_resource
    # * submit -> submit -> submitting
    # * wait_until_execution_starts -> wait_until_execution_starts -> waiting_for_execution
    # * wait_for_execution_end -> wait_for_execution_end -> executing
    # * transfer_data_from -> transfer_data_from_resource -> transferring_data_from_resource
    # * store_result_in_backend -> store_result_in_backend -> importing_into_database
    # * cleanup -> cleanup -> cleaning_up_workspace
    # * success -> success -> completed_successfully
    # * error -> error -> completed_error
    # * <from PUT status> -> cancelling -> cancelling
    # * cancel -> cancelled *RENOMBRADO* -> cancelled
    # * ... -> ... -> paused
    status = Column(String(80))
    log = Column(Text)
    outputs = Column(JSON)
    # To store (and be able to resume) "task" (in definitions.py) and job_context (needed to follow)
    execution_state = Column(JSONB)
    deleted = Column(Boolean, default=False)  # True: soft-deleted


# SUBSYSTEM STATUS

prefix = "sa_status_"


class StatusEnum(Enum):
    working = 'This system works.'
    failing = 'This system does not work.'
    undefined = 'System status is unknown.'

    def __str__(self):
        return json.dumps(self, default=lambda x: x.value)


class StatusChecker(ORMBase):
    __tablename__ = f"{prefix}checkers"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    name = Column(String(32), nullable=False)
    type = Column(String(32), nullable=False)   # pg, redis, ssh, compute-resources
    url = Column(String(80))
    resource_id = Column(Integer, ForeignKey(ComputeResource.id), nullable=True, primary_key=False)
    resource = relationship(ComputeResource)
    status = Column(ENUM(StatusEnum))
    details = Column(JSON)
    last_check = Column(DateTime, default=datetime.utcnow())

    __table_args__ = (
        UniqueConstraint(type, url, name=__tablename__ + '_c1'),
        UniqueConstraint(type, resource_id, name=__tablename__ + '_c2'),
        CheckConstraint('COALESCE(url , resource_id::text) IS NOT NULL', name=__tablename__ + '_c3'),
    )

    def check_status(self):
        from urllib.parse import urlparse
        u = urlparse(self.url)
        if self.type == 'compute-resource':
            from ..jobs.ssh_resource import JobExecutorWithSSH
            je = JobExecutorWithSSH('0', create_local_workspace=False, initialize_loop=False)
            je.set_resource(json.loads(generate_json(self.resource)))
            self.details = _status = je.check()
        elif self.type == 'redis' or self.type == 'redis-cli':
            # cmd = 'redis-cli -h %s -p %s ping' % (u.hostname, u.port)
            from redis import Redis, ConnectionError
            try:
                self.details = Redis(u.hostname, u.port).ping()
                _status = True
            except ConnectionError:
                _status = False
        else:   # cmd with subprocess
            if self.type == 'pg' or u.scheme == 'postgres' or u.scheme == 'postgresql':
                cmd = 'pg_isready'
                cmd += ' -h %s' % u.hostname
                cmd += ' -p %s' % u.port if u.port else ''
                cmd += ' -U %s' % u.username if u.username else ''
                cmd += ' -d %s' % u.path[1:] if u.path[1:] else ''
            else:   # elif self.type == 'curl' or u.scheme == 'http' or u.scheme == 'https':
                from ..services import secure_url
                cmd = 'curl --fail -s %s || exit 1' % secure_url(self.url)
            import subprocess
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            out, err = process.communicate()
            self.details = {'stdout': out.decode('utf-8'), 'stderr': err.decode('utf-8')}
            _status = not process.returncode
        self.status = StatusEnum.undefined if _status is None else StatusEnum.working if _status else StatusEnum.failing
        self.last_check = datetime.utcnow()
        return _status

