import asyncio
import os
import psutil
from shutil import rmtree

import asyncssh

from .ssh_process_adaptors import SSHProcessAdaptor
from .. import get_global_configuration_variable
from ..jobs import JobExecutorAtResource


class CustomSSHClient(asyncssh.SSHClient):

    def validate_host_public_key(self, host, addr, port, key):
        print(f'host: {host}, addr: {addr}, port: {port}, key: {key}')
        return True


class RemoteSSHClient:
    SSH_OPTIONS = "-o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no"
    RSYNC_OPTIONS = f"-e 'ssh {SSH_OPTIONS}'"
    """Client to interact with a remote host via SSH & SCP."""

    def __init__(self, host, data_host, port, data_port, username, known_hosts_filepath, remote_workspace, local_workspace, logs_dict):
        super()
        self.host = host
        self.data_host = data_host
        self.port = port
        self.data_port = data_port
        self.username = username
        self.known_hosts_filepath = known_hosts_filepath
        self.remote_workspace = remote_workspace
        self.local_workspace = local_workspace
        self.logs_dict = logs_dict
        self.client = None
        self.sftp = None
        self.conn = None

    async def connect(self):
        """
        Open connection to remote host.
        @return: RemoteSSHClient connected instance
        """
        if self.conn is None:
            try:
                kwargs = {  # I need known_hosts and the arguments in the constructor of RemoteSSHClient
                    'known_hosts': self.known_hosts_filepath,
                }
                self.conn, self.client = await asyncssh.create_connection(CustomSSHClient,
                                                                          self.host, self.port, username=self.username)
                self.sftp = await self.conn.start_sftp_client()
                # NOTE: >>> the local workspace is created by JobExecutorAtResource.__init__ <<<
                # if create_local_workspace:
                #     if not os.path.exists(self.local_workspace):  # create folder for saving local jobs exit status
                #         os.mkdir(self.local_workspace)
            except Exception as error:
                print(f'Connection failed: \
                   did you remember to create a known_host file? {error}')
                raise error
        return self

    def disconnect(self):
        """Close ssh connection."""
        if self.sftp:
            self.sftp.exit()
        if self.conn:
            self.conn.close()

    async def run_client(self, script_file, script_params):
        """
        Execute remote client script
        @param script_file: Path to the script file
        @param script_params: Parameters of the script
        @return: pid: PID of the executed script process
        """
        cmd = (f"ssh {self.SSH_OPTIONS} -p {self.port} {self.username}@{self.host} 'cd {self.remote_workspace} && export PATH=.:$PATH && chmod +x {script_file} &> /dev/null " +
               f"; (nohup {script_file} {script_params} " +
               f">/{self.remote_workspace}/{os.path.basename(self.logs_dict['submit_stdout'])} " +
               f"</dev/null 2>/{self.remote_workspace}/{os.path.basename(self.logs_dict['submit_stderr'])}" +
               f"& echo $!; wait $!; echo $? >> {self.remote_workspace}/$!.exit_status)'")

        print(repr(cmd))
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().rstrip()
        print(f"PID: {pid}")
        return pid

    async def write_submit_logs(self):
        remote_stdout_file = await self.sftp.open((f"{self.remote_workspace}" +
                                                   f"/{os.path.basename(self.logs_dict['submit_stdout'])}"))
        stdout_content = await remote_stdout_file.read()
        await remote_stdout_file.close()
        with open(self.logs_dict['submit_stdout'], "a") as file:
            file.write(stdout_content)

        remote_stderr_file = await self.sftp.open((f"{self.remote_workspace}" +
                                                   f"/{os.path.basename(self.logs_dict['submit_stderr'])}"))
        stderr_content = await remote_stderr_file.read()
        await remote_stderr_file.close()
        with open(self.logs_dict['submit_stderr'], "a") as file:
            file.write(stderr_content)

    async def remote_command_status(self, pid):
        """
        Get Job status
        @param pid: PID of the process to check status
        @return: String defining satus of the pid. "running", "ok" and "" for error.
        """
        if await self.exists_remotely(f"{self.remote_workspace}/{pid}.exit_status"):
            cmd = f"ssh {self.SSH_OPTIONS} -p {self.port} {self.username}@{self.host} 'cat {self.remote_workspace}/{pid}.exit_status'"
            popen_pipe = os.popen(cmd)
            exit_status = popen_pipe.readline().strip()
            if exit_status.strip() == "0":
                exit_status = "ok"
            else:
                print(
                    f"Error executing remote job with pid: {pid}. Exit status = {exit_status}; Host = {self.host}")
                exit_status = ""  # This means error
        else:
            exit_status = "running"
        print(f"Job Status: {exit_status}")

        return exit_status

    def local_command_status(self, pid):
        """
        Get Job status
        @param pid: PID of the process to check status
        @return: String defining satus of the pid. "running", "ok" and "" for error.
        """
        if os.path.exists(f"{self.local_workspace}/{pid}.exit_status"):
            cmd = f"cat {self.local_workspace}/{pid}.exit_status"
            popen_pipe = os.popen(cmd)
            exit_status = popen_pipe.readline().strip()
            if exit_status.strip() == "0":
                exit_status = "ok"
            else:
                print(f"Error executing local job with pid: {pid}. Exit status = {exit_status}")
                exit_status = ""  # This means error
        else:
            exit_status = "running"

        print(f"Job Status: {exit_status}")
        return exit_status

    def kill_process(self, pid):
        """
        @param pid: PID of the process to kill
        """
        if pid is not None and pid != "":
            last_job_remotely = not psutil.pid_exists(int(pid))
            if last_job_remotely:
                print(f"Cancel command: ssh {self.SSH_OPTIONS} -p {self.port} {self.username}@{self.host} 'kill -9 -- -$(ps -p {pid} -o pgid=)'")
                os.system(f"ssh {self.SSH_OPTIONS} -p {self.port} {self.username}@{self.host} 'kill -9 -- -$(ps -p {pid} -o pgid=)'")
            else:
                print(
                    f"Cancel command: kill -9 -- -$(ps -p {pid} -o pgid=)")
                os.system(f"kill -9 -- -$(ps -p {pid} -o pgid=)")
        else:
            print("No command has been executed")

    def upload_file(self, local_path, remote_path):
        """
        Upload file to remote_path.

        @param local_path: path to local file.
        @param remote_path: the remote path to file relative to self.remote_workspace.
        @return: PID of the process uploading the file.
        """
        remote_path = os.path.join(self.remote_workspace, remote_path)
        cmd = (f"(nohup scp -P {self.port} {self.SSH_OPTIONS} {local_path} {self.username}@{self.host}:{remote_path} " +
                   f">>{self.logs_dict['download_stdout']} </dev/null 2>>{self.logs_dict['download_stderr']} " +
                   f"& echo $!; wait $!; echo $? >> {self.local_workspace}/$!.exit_status)")
        print(cmd)
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().rstrip()
        print(f"PID: {pid}")
        return pid

    async def download_file(self, remote_file, local_file):
        """Download file from remote host.
        @param remote_file: path to remote file relative to self.remote_workspace.
        @param local_file: local path name of downloaded file.
        @return: PID of the process downloading the file.
        """
        pid = None
        if self.sftp is not None:
            remote_file = os.path.join(self.remote_workspace, remote_file)
            if await self.sftp.isfile(remote_file):
                cmd = (f"(nohup scp -P {self.port} {self.SSH_OPTIONS} {self.username}@{self.host}:{remote_file} {local_file} " +
                       f">>{self.logs_dict['download_stdout']} </dev/null 2>>{self.logs_dict['download_stderr']} " +
                       f"& echo $!; wait $!; echo $? >> {self.local_workspace}/$!.exit_status)")

                print(cmd)
                popen_pipe = os.popen(cmd)
                pid = popen_pipe.readline().rstrip()
                print(f"PID: {pid}")
            else:
                print("The remote file doesn't exist")
        else:
            print("SSH connection not created")
        return pid

    async def move_file(self, filepath, new_filepath):
        """
        Move remote file in the host
        @param filepath: Current path of the remote path
        @param new_filepath: New path of the remote path
        """
        if self.sftp is not None:
            await self.sftp.rename(os.path.join(self.remote_workspace, filepath),
                                   os.path.join(self.remote_workspace, new_filepath))
        else:
            print("SSH connection not created")

    async def remove_file(self, file):
        """
        Remove remote path
        @param file: Path to the remote file
        """
        if self.sftp is not None:
            await self.sftp.remove(os.path.join(self.remote_workspace, file))
        else:
            print("SSH connection not created")

    async def make_directory(self, dir_name):
        """
        Make a remote directory
        @param dir_name: Name of the remote directory relative to self.remote_workspace.
        If an absolute path is given it is not relative to self.remote_workspace.
        """
        if self.sftp is not None:
            if not await self.sftp.exists(os.path.join(self.remote_workspace, dir_name)):
                path = self.remote_workspace if self.remote_workspace == dir_name else os.path.join(self.remote_workspace, dir_name)
                await self.sftp.mkdir(path)
            else:
                print("folder already exists")
        else:
            print("SSH connection not created")

    async def directory_exists(self, dir_name):
        """
        Check remote directory exists
        @param dir_name: Name of the remote directory relative to self.remote_workspace.
        If an absolute path is given it is not relative to self.remote_workspace.
        """
        return await self.sftp.exists(os.path.join(self.remote_workspace, dir_name))

    async def remove_directory(self, dir_path):
        """
        Remove remote directory
        @param dir_path: Remote directory path relative to self.remote_workspace
        """
        if self.sftp is not None:
            if await self.sftp.isdir(os.path.join(self.remote_workspace, dir_path)):
                await self.sftp.rmtree(os.path.join(self.remote_workspace, dir_path))
            else:
                print(
                    (f"The path {os.path.join(self.remote_workspace, dir_path)} " +
                     "doesn't correspond to a remote directory."))
        else:
            print("SSH connection not created")

    async def exists_remotely(self, name):
        """
        Exists remote path
        @param name: Remote path
        @return: Boolean value indicating if path exists in host
        """
        print(f"Exists remotely check: {os.path.join(self.remote_workspace, name)}")
        exist = await self.sftp.exists(os.path.join(self.remote_workspace, name))
        return exist

    async def same_size(self, local_name, remote_name):
        """
        Compares the local path size with the remote path size.
        It ONLY works with files and empty directories
        @param local_name: Local Path
        @param remote_name: Remote Path
        @return: Boolean value indicating that they have the same size
        """
        same_size = os.path.getsize(local_name) == await self.sftp.getsize(
            os.path.join(self.remote_workspace, remote_name))
        if not same_size:
            print(f"The files {local_name} and {remote_name} don't have the same size")
        return same_size


class JobExecutorWithSSH(JobExecutorAtResource):

    def __init__(self, identity_job_id, create_local_workspace=True, remote_workspace=None, initialize_loop=True):
        super().__init__(identity_job_id, create_local_workspace)
        self.host = None
        self.data_host = None
        self.port = 22
        self.data_port = 22
        self.username = None
        self.known_hosts_filepath = None
        self.remote_workspace = os.path.join(remote_workspace, identity_job_id) if remote_workspace else os.path.join(get_global_configuration_variable("SSH_JOBS_DEFAULT_REMOTE_WORKSPACE"),
                                             identity_job_id)
        self.remote_client = None
        self.loop = asyncio.get_event_loop() if initialize_loop else None

    # RESOURCE
    def set_resource(self, resource_params):
        self.host = resource_params["jm_location"]['host']
        self.port = int(resource_params["jm_location"].get('port', "22"))
        self.data_host = resource_params["jm_location"].get('data_host', self.host)
        self.data_port = resource_params["jm_location"].get('data_port', self.port)
        self.username = resource_params["jm_credentials"]['username']
        self.known_hosts_filepath = resource_params["jm_credentials"]['known_hosts_filepath']

    def check(self):
        not_accessible = os.system(f"ssh -p {self.port} -o StrictHostKeyChecking=no -o " +
                                   f"UserKnownHostsFile={self.known_hosts_filepath} {self.username}@{self.host} " +
                                   f"'echo'")
        print(not_accessible)
        if not_accessible:  # if not_accessible is 0 the server is accessible else it is not
            return False
        else:
            return True

    def connect(self):
        self.remote_client = RemoteSSHClient(self.host, self.data_host, self.port, self.data_port, self.username,
                                             self.known_hosts_filepath, self.remote_workspace,
                                             self.local_workspace, self.log_filenames_dict)
        self.loop.run_until_complete(self.remote_client.connect())
        return self.remote_client

    def disconnect(self):
        self.remote_client.disconnect()
        self.loop.close()

    # JOB EXECUTION
    def create_job_workspace(self):
        # the name is the job_id
        self.loop.run_until_complete(self.remote_client.make_directory(self.remote_workspace))

    def job_workspace_exists(self):
        cmd = f"ssh {self.remote_client.SSH_OPTIONS} -p {self.port} -G {self.host}" + " | awk 'FNR == 2 {print $2}'"
        popen_pipe = os.popen(cmd)
        host_ip = popen_pipe.readline().rstrip()
        if host_ip == "127.0.0.1":
            return os.path.isdir(self.remote_workspace)
        _ = self.loop.run_until_complete(self.remote_client.directory_exists(""))
        return _

    def remove_job_workspace(self):
        cmd = f"ssh {self.remote_client.SSH_OPTIONS} -p {self.port} -G {self.host}" + " | awk 'FNR == 2 {print $2}'"
        popen_pipe = os.popen(cmd)
        host_ip = popen_pipe.readline().rstrip()
        if host_ip == "127.0.0.1":
            return rmtree(self.remote_workspace)
        else:
            # cmd = ["ssh", f"{SSH_OPTIONS}",
            #        "-p", f"{self.port}", f"{self.username}@{self.host}",
            #        f"rm -fr {self.remote_workspace}"]
            # proc = subprocess.run(cmd, capture_output=True, text=True)
            # print(f"STDOUT: {proc.stdout}")
            # print(f"STDERR: {proc.stderr}")
            self.loop.run_until_complete(self.remote_client.remove_directory(""))

    def exists(self, job_context):
        i = job_context["state_dict"]["idx"]
        if job_context["state_dict"]["state"] == "upload":
            local_path = self.get_upload_files_list(job_context)[i]["file"]
            remote_path = self.get_upload_files_list(job_context)[i]["remote_name"]
        else:
            filename = self.get_download_files_list(job_context)[i]["file"]
            local_path = os.path.join(self.local_workspace, filename)
            remote_path = self.get_download_files_list(job_context)[i]["remote_name"]

        if os.path.exists(local_path):
            check = self.loop.run_until_complete(self.remote_client.exists_remotely(remote_path))
            if check:
                # the same_size method works ONLY with files and empty directories
                check &= self.loop.run_until_complete(
                    self.remote_client.same_size(local_path, remote_path))
            else:
                print(f"File {remote_path} not found in the remote system")
            return check
        else:
            print(f"File {local_path} not found in your local system")
            return None

    def upload_file(self, job_context):
        i = job_context["state_dict"]["idx"]
        local_path = self.get_upload_files_list(job_context)[i]["file"]
        remote_path = self.get_upload_files_list(job_context)[i]["remote_name"]
        return self.remote_client.upload_file(local_path, remote_path)

    def download_file(self, job_context):
        i = job_context["state_dict"]["idx"]
        remote_path = self.get_download_files_list(job_context)[i]["remote_name"]
        filename = self.get_download_files_list(job_context)[i]["file"]
        local_path = os.path.join(self.local_workspace, filename)
        return self.loop.run_until_complete(self.remote_client.download_file(remote_path, local_path))

    def submit(self, process):
        params = process["inputs"]["adapted_parameters"]
        return self.loop.run_until_complete(
            self.remote_client.run_client(params[SSHProcessAdaptor.SCRIPT_KEY][0]["remote_name"],
                                          params[SSHProcessAdaptor.SCRIPT_PARAMS_KEY]))

    def step_status(self, job_context):
        pid = job_context["pid"]
        previous_state = job_context["state_dict"]["state"]

        if pid is None:
            status = "none"
        elif pid == "":
            status = ""  # error
        elif previous_state == "submit":
            status = self.loop.run_until_complete(self.remote_client.remote_command_status(pid))
        else:
            status = self.local_job_status(pid)

        if status == "ok" or status == "":
            self.write_remote_logs(job_context["state_dict"])

        return status

    def write_remote_logs(self, state_dict):
        if state_dict["state"] == "submit":
            self.loop.run_until_complete(self.remote_client.write_submit_logs())

    def cancel_job(self, native_id):
        self.remote_client.kill_process(native_id)

    def get_upload_files_list(self, job_context):
        # TODO : check if script files appear duplicated in the resulting list
        return list(job_context["process"]["inputs"]["data"]) + \
               list(job_context["process"]["inputs"]["adapted_parameters"][SSHProcessAdaptor.SCRIPT_FILES_KEY]) + \
               list(job_context["process"]["inputs"]["adapted_parameters"].get("scripts", []))

    def get_download_files_list(self, job_context):
        return job_context["results"]

    def check_resource(self):
        """
        Can connect
        CPU and GPU use
        Available storage space

        :return:
        """
        # TODO:
        pass
