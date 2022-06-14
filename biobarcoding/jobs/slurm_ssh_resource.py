import os

import psutil

from .ssh_resource import RemoteSSHClient, JobExecutorWithSSH
from .. import get_global_configuration_variable


class RemoteSlurmClient(RemoteSSHClient):

    JOB_STATES_DICT = {
        'BOOT_FAIL': "",# "" means error
        'CANCELLED': "",
        'COMPLETED': "ok",
        'CONFIGURING': "wait_until_start",
        'COMPLETING': "running",
        'DEADLINE': "",
        'FAILED': "",
        'NODE_FAIL': "",
        'OUT_OF_MEMORY': "",
        'PENDING': "wait_until_start",
        'PREEMPTED': "running",#TODO ?
        'RUNNING': "running",
        'RESV_DEL_HOLD': "running",#TODO ?
        'REQUEUE_FED': "running",
        'REQUEUE_HOLD': "running",
        'REQUEUED': "running",
        'RESIZING': "running",
        'REVOKED': "running",
        'SIGNALING': "running",
        'SPECIAL_EXIT': "running",#TODO ?
        'STAGE_OUT': "running",
        'STOPPED': "running",#TODO ?
        "SUSPENDED": "running",#TODO ?
        "TIMEOUT": ""
    }

    async def run_client(self, script_file, script_params):
        """
        Execute remote client script
        @param script_file: Path to the script file
        @param script_params: Parameters of the script and hpc (ncpus, ngpus, time..)
        @return: pid: Job ID of the executed script process
        """
        hpc_parameters = script_params["hpc_parameters"]
        ntasks = f"--ntasks={hpc_parameters['ntasks']}" if 'ntasks' in hpc_parameters else ""
        cpus_per_task = f"--cpus-per-task={hpc_parameters['cpus_per_task']}" if 'cpus_per_task' in hpc_parameters else ""
        cpus_per_gpu = f"--cpus-per-gpu={hpc_parameters['cpus_per_gpu']}" if 'cpus_per_gpu' in hpc_parameters else ""
        gpus = f"--gpus={hpc_parameters['gpus']}" if 'gpus' in hpc_parameters else ""
        '''mem_per_gpu = f"--mem-per-gpu={hpc_parameters['mem_per_gpu']}" if 'mem_per_gpu' in hpc_parameters else ""
        mem_per_cpu = f"--mem-per-cpu={hpc_parameters['mem_per_cpu']}" if 'mem_per_cpu' in hpc_parameters else ""'''
        priority = f"--priority={hpc_parameters['priority']}" if 'priority' in hpc_parameters else ""
        #threads_per_core = f"--threads-per-core={hpc_parameters['threads_per_core']}" if 'threads_per_core' in hpc_parameters else ""
        time = f"--time={hpc_parameters['time']}" if 'time' in hpc_parameters else ""#TODO hablar lo del PENDING
        '''time_min = f"--time-min={hpc_parameters['time_min']}" if 'time_min' in hpc_parameters else ""
        tmp = f"--tmp={hpc_parameters['tmp']}" if 'tmp' in hpc_parameters else ""'''

        process_parameters = script_params['process_parameters']

        cmd = (
                f"ssh {self.SSH_OPTIONS} {self.username}@{self.host} \"cd {self.remote_workspace} " +
                f"&& sbatch {ntasks} {cpus_per_task} {cpus_per_gpu} {gpus} {priority} " +
                f"{time} --export=ALL,{process_parameters} --chdir={self.remote_workspace} " +
                f"--open-mode=truncate --error={self.remote_workspace}/{os.path.basename(self.logs_dict['submit_stderr'])} " +
                f"--output={self.remote_workspace}/{os.path.basename(self.logs_dict['submit_stdout'])} " +
                f"{script_file}\"" + " | awk 'NF{print $NF; exit}'")

        print(repr(cmd))
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().strip()
        print(f"PID: {pid}")
        return pid

    async def remote_command_status(self, native_job_id):
        """
        Get Job status
        @param native_job_id: Job ID of the process to check status
        @return: String defining satus of the pid. "running", "ok" and "" for error.
        """
        cmd = f"ssh {self.SSH_OPTIONS} {self.username}@{self.host} \"scontrol show job {native_job_id}\" | grep 'JobState' | sed 's/=/ /' | awk '{{print $2}}'"
        popen_pipe = os.popen(cmd)
        job_state = popen_pipe.readline().strip()
        if job_state != '':
            exit_status = self.JOB_STATES_DICT[job_state]
            print(f"Job State: {job_state}")
            print(f"Exit Status: {exit_status}")
            return exit_status
        else:
            return "running"

    def kill_process(self, pid):
        """
        @param pid: Job ID to kill
        """
        if pid is not None and pid != "":
            last_job_remotely = psutil.pid_exists(int(pid))
            if last_job_remotely:
                os.system(f"ssh {self.username}@{self.host} 'scancel {pid}'")
            else:
                os.system(f"kill -9 -- -$(ps -p {pid} -o pgid=)")
        else:
            print("No command has been executed")


class JobExecutorWithSlurm(JobExecutorWithSSH):

    def connect(self):
        self.remote_client = RemoteSlurmClient(self.host, self.data_host, self.port, self.data_port, self.username,
                                               self.known_hosts_filepath, self.remote_workspace,
                                               self.local_workspace, self.log_filenames_dict)
        self.loop.run_until_complete(self.remote_client.connect())
        return self.remote_client
