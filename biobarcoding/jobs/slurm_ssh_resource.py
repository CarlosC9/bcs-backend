import os
from .ssh_resource import JobExecutorWithSSH, RemoteSSHClient
from .. import get_global_configuration_variable


class RemoteSlurmClient(RemoteSSHClient):

    JOB_STATES_DICT = {
        'BOOT_FAIL': "",# "" means error
        'CANCELLED': "",
        'COMPLETED': "ok",
        'CONFIGURING': "running",
        'COMPLETING': "running",
        'DEADLINE': "",
        'FAILED': "",
        'NODE_FAIL': "",
        'OUT_OF_MEMORY': "",
        'PENDING': "running",
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
        resources_params = script_params["resources_params"]
        cpus_per_task = f"--cpus-per-task={resources_params['cpus_per_task']}" if 'cpus_per_task' in resources_params else ""
        cpus_per_gpu = f"--cpus-per-gpu={resources_params['cpus_per_gpu']}" if 'cpus_per_gpu' in resources_params else ""
        gpus = f"--gpus={resources_params['gpus']}" if 'gpus' in resources_params else ""
        mem_per_gpu = f"--mem-per-gpu={resources_params['mem_per_gpu']}" if 'mem_per_gpu' in resources_params else ""
        mem_per_cpu = f"--mem-per-cpu={resources_params['mem_per_cpu']}" if 'mem_per_cpu' in resources_params else ""
        priority = f"--priority={resources_params['priority']}" if 'priority' in resources_params else ""
        threads_per_core = f"--threads-per-core={resources_params['threads_per_core']}" if 'threads_per_core' in resources_params else ""
        time = f"--time={resources_params['time']}" if 'time' in resources_params else ""#TODO hablar lo del PENDING
        time_min = f"--time-min={resources_params['time_min']}" if 'time_min' in resources_params else ""
        tmp = f"--tmp={resources_params['tmp']}" if 'tmp' in resources_params else ""
        process_params = ""
        for key, value in script_params['process_params']:
            process_params += f"{key}={value},"
        process_params = process_params[:-1] #remove last comma
        cmd = (
                f"ssh {self.SSH_OPTIONS} {self.username}@{self.host} 'cd {self.remote_workspace} " +
                f"&& sbatch {cpus_per_task} {cpus_per_gpu} {gpus} {mem_per_gpu} {mem_per_cpu} {priority} " +
                f"{threads_per_core} {time} {time_min} {tmp} --export=ALL,{process_params} --chdir={self.remote_workspace} " +
                f"--open-mode=append --error={self.remote_workspace}/{os.path.basename(self.logs_dict['submit_stderr'])} " +
                f"--output={self.remote_workspace}/{os.path.basename(self.logs_dict['submit_stdout'])} " +
                f"{script_file}" + " | awk 'NF{print $NF; exit}'")

        print(repr(cmd))
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().strip()
        self.last_job_remotely = True
        print(f"PID: {pid}")
        return pid

    async def remote_command_status(self, pid):
        """
        Get Job status
        @param pid: Job ID of the process to check status
        @return: String defining satus of the pid. "running", "ok" and "" for error.
        """
        cmd = f"ssh {self.SSH_OPTIONS} {self.username}@{self.host} 'sacct -n -X --jobs={pid} --format=state'"
        popen_pipe = os.popen(cmd)
        job_state = popen_pipe.readline().strip()
        exit_status = self.JOB_STATES_DICT[job_state]
        print(f"Job State: {job_state}")
        print(f"Exit Status: {exit_status}")
        return exit_status

    def kill_process(self, pid):
        """
        @param pid: Job ID to kill
        """
        if pid is not None and pid != "":
            if self.last_job_remotely:
                os.system(f"ssh {self.username}@{self.host} 'scancel {pid}'")
            else:
                os.system(f"kill -9 {pid}")
        else:
            print("No command has been executed")


class JobExecutorWithSlurm(JobExecutorWithSSH):

    def connect(self):
        self.remote_client = RemoteSlurmClient(self.host, self.username,
                                               self.known_hosts_filepath, self.remote_workspace,
                                               self.local_workspace, self.log_filenames_dict)
        self.loop.run_until_complete(self.remote_client.connect())
        return self.remote_client

    def submit(self, process):
        params = process["inputs"]["parameters"]
        return self.loop.run_until_complete(
            self.remote_client.run_client(params["scripts"][0]["remote_name"],
                                          {"process_params": params["script_params"],
                                           "resources_params": params["hpc_params"]}))
