import psutil

from biobarcoding.jobs import JobExecutorAtResource
import os
import asyncio
import asyncssh

#TODO revise this
ROOT_DIR = "/home/dreyes"


class RemoteSSHClient:
    """Client to interact with a remote host via SSH & SCP."""

    def __init__(self, host, username, known_hosts_filepath, remote_path):
        super()
        self.host = host
        self.username = username
        self.known_hosts_filepath = known_hosts_filepath
        self.remote_path = remote_path
        self.job_status_dir = os.path.join("/tmp", username + "_jobs_status")
        self.client = None
        self.sftp = None
        self.conn = None
        self.last_job_remotely = None

    async def connect(self):
        """
        Open connection to remote host.
        @return: RemoteSSHClient connected instance
        """
        if self.conn is None:
            try:
                kwargs = {  # necesito known_hosts y los argumentos del constructor de RemoteSSHClient
                    'known_hosts': self.known_hosts_filepath,
                }
                self.conn, self.client = await asyncssh.create_connection(asyncssh.SSHClient, self.host, **kwargs)
                self.sftp = await self.conn.start_sftp_client()
                if not os.path.isdir(self.job_status_dir): #create folder for saving local jobs exit status
                    os.mkdir(self.job_status_dir)
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
        job_id = os.path.basename(os.path.normpath(self.remote_path))
        cmd = f"ssh {self.username}@{self.host} 'cd {self.remote_path} && chmod +x {script_file} && (nohup ./{script_file} {script_params}  >/{self.remote_path}/bcs.stdout.log </dev/null 2>/{self.remote_path}/bcs.err.log & echo $!; wait $!; echo $? >> {self.remote_path}/$!.exit_status)'"

        print(cmd)
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().rstrip()
        self.last_job_remotely = True
        print(f"PID: {pid}")
        return pid

    async def command_status(self, pid):
        """
        Get Job status
        @param pid: PID of the process to check status
        @return: String defining satus of the pid. "running", "ok" and "" for error.
        """
        exit_status = ""#This means error
        if pid is not None and pid != "":
            if self.last_job_remotely:
                if await self.exists_remotely(f"{self.remote_path}/{pid}.exit_status"):
                    cmd = f"ssh {self.username}@{self.host} 'cat {self.remote_path}/{pid}.exit_status'"
                    popen_pipe = os.popen(cmd)
                    exit_status = popen_pipe.readline().strip()
                    if exit_status.strip() == "0":
                        exit_status = "ok"
                    else:
                        print(f"Error executing remote job with pid: {pid}. Exit status = {exit_status}; Host = {self.host}")
                        exit_status = ""#This means error
                else:
                    exit_status = "running"
            else:
                if os.path.isfile(f"{self.job_status_dir}/{pid}.exit_status"):
                    cmd = f"cat {self.job_status_dir}/{pid}.exit_status"
                    popen_pipe = os.popen(cmd)
                    exit_status = popen_pipe.readline().strip()
                    if exit_status.strip() == "0":
                        exit_status = "ok"
                    else:
                        print(f"Error executing local job with pid: {pid}. Exit status = {exit_status}")
                        exit_status = ""#This means error
                else:
                    exit_status = "running"
            print(f"Exit Status: {exit_status}")
        else:
            print(f"The PID can't be None")
        return exit_status

    def kill_process(self, pid):
        """
        @param pid: PID of the process to kill
        """
        if pid is not None and pid != "":
            if self.last_job_remotely:
                os.system(f"ssh {self.username}@{self.host} 'kill -9 {pid}'")
            else:
                os.system(f"kill -9 {pid}")
        else:
            print("No command has been executed")

    def upload_directory(self, local_path, remote_path):
            """
            Upload directory to remote_path.

            @param local_path: path to local folder.
            @param remote_path: path to remote folder relative to self.remote_path.
            @return: PID of the process uploading the directory
            """
            remote_path = os.path.join(self.remote_path, remote_path, "")#ensure it always ends with / (the separator)
            remote_dir = os.path.join(self.remote_path, os.path.split(remote_path)[0])
            if remote_dir != "" and False:
                create_remote_dir_cmd = f"--rsync-path='mkdir -p {remote_dir} & rsync'"
            else:
                create_remote_dir_cmd = ""
            #-a for folders
            cmd = f"(nohup bash -c \"rsync -a {create_remote_dir_cmd} {local_path} {self.username}@{self.host}:{remote_path}\" >/tmp/mtest2 </dev/null 2>/tmp/mtest2.err & echo $!; wait $!; echo $? >> {self.job_status_dir}/$!.exit_status)"
            print(cmd)
            popen_pipe = os.popen(cmd)
            self.last_job_remotely = False
            pid = popen_pipe.readline().rstrip()
            print(f"PID: {pid}")
            return pid

    async def download_directory(self, remote_dir, local_dir):
        """
        Download directory from remote host.
        @param remote_dir: path to remote folder relative to self.remote_path.
        @param local_dir: local path name of downloaded folder.
        @return: PID of the process downloading the directory.
        """
        pid = None
        if self.sftp is not None:
            remote_dir = os.path.normpath(os.path.join(self.remote_path, remote_dir)) #ensure it doesn't finish with / (separator)
            if self.sftp.isdir(os.path.join(self.remote_path, remote_dir)):
                local_dir = os.path.join(local_dir, "")#ensure that it ends with / (separator)
                cmd = f"(nohup \"scp -r {self.username}@{self.host}:{remote_dir} {local_dir}\" (>/tmp/mtest2 </dev/null 2>/tmp/mtest2.err & echo $!; wait $!; echo $? >> {self.job_status_dir}/$!.exit_status))"
                print(cmd)
                popen_pipe = os.popen(cmd)
                self.last_job_remotely = False
                pid = popen_pipe.readline().rstrip()
                print(f"PID: {pid}")
            else:
                print("The remote directory doesn't exist")
        else:
            print("SSH connection not created")
        return pid

    def upload_file(self, local_path, remote_path):
            """
            Upload file to remote_path.

            @param local_path: path to local file.
            @param remote_path: the remote path to file relative to self.remote_path.
            @return: PID of the process uploading the file.
            """
            remote_path = os.path.join(self.remote_path, remote_path)
            remote_dir = os.path.join(self.remote_path, os.path.split(remote_path)[0])
            if remote_dir != "" and False:
                create_remote_dir_cmd = f"--rsync-path='mkdir -p {remote_dir} & rsync'"
            else:
                create_remote_dir_cmd = ""
            cmd = f"(nohup bash -c \"rsync {create_remote_dir_cmd} {local_path} {self.username}@{self.host}:{remote_path}\" >/tmp/mtest2 </dev/null 2>/tmp/mtest2.err & echo $!; wait $!; echo $? >> {self.job_status_dir}/$!.exit_status)"
            print(cmd)
            popen_pipe = os.popen(cmd)
            self.last_job_remotely = False
            pid = popen_pipe.readline().rstrip()
            print(f"PID: {pid}")
            return pid

    async def download_file(self, remote_file, local_file):
        """Download file from remote host.
        @param remote_file: path to remote file relative to self.remote_path.
        @param local_file: local path name of downloaded file.
        @return: PID of the process downloading the file.
        """
        pid = None
        if self.sftp is not None:
            remote_file = os.path.join(self.remote_path, remote_file)
            if self.sftp.isfile(remote_file):
                cmd = f"(nohup \"scp {self.username}@{self.host}:{remote_file} {local_file}\" (>/tmp/mtest2 </dev/null 2>/tmp/mtest2.err & echo $!; wait $!; echo $? >> {self.job_status_dir}/$!.exit_status))"
                print(cmd)
                popen_pipe = os.popen(cmd)
                self.last_job_remotely = False
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
            await self.sftp.rename(os.path.join(self.remote_path, filepath),
                                   os.path.join(self.remote_path, new_filepath))
        else:
            print("SSH connection not created")

    async def remove_file(self, file):
        """
        Remove remote path
        @param file: Path to the remote file
        """
        if self.sftp is not None:
            await self.sftp.remove(os.path.join(self.remote_path, file))
        else:
            print("SSH connection not created")

    async def make_directory(self, dir_name):
        """
        Make a remote directory
        @param dir_name: Name of the remote directory relative to self.remote_path.
        If an absolute path is given it is not relative to self.remote_path.
        """
        if self.sftp is not None:
            if not await self.sftp.isdir(os.path.join(self.remote_path, dir_name)):
                # TODO: Se le puede asignar permisos, por ejemplo al working directory
                await self.sftp.mkdir(os.path.join(self.remote_path, dir_name))
            else:
                print("folder already exists")
        else:
            print("SSH connection not created")

    async def remove_directory(self, dir_path):
        """
        Remove remote directory
        @param dir_path: Remote directory path relative to self.remote_path
        """
        if self.sftp is not None:
            if self.sftp.isdir(os.path.join(self.remote_path, dir_path)):
                await self.sftp.rmtree(os.path.join(self.remote_path, dir_path))
            else:
                print(f"The path {os.path.join(self.remote_path, dir_path)} doesn't correspond to a remote directory.")
        else:
            print("SSH connection not created")

    async def exists_remotely(self, name):
        """
        Exists remote path
        @param name: Remote path
        @return: Boolean value indicating if path exists in host
        """
        print(f"Exists remotely: {os.path.join(self.remote_path, name)}")
        return await self.sftp.exists(os.path.join(self.remote_path, name))

    async def same_size(self, local_name, remote_name):
        """
        Compares the local path size with the remote path size.
        It ONLY works with files and empty directories
        @param local_name: Local Path
        @param remote_name: Remote Path
        @return: Boolean value indicating that they have the same size
        """
        return os.path.getsize(local_name) == await self.sftp.getsize(os.path.join(self.remote_path, remote_name))


class JobExecutorWithSSH(JobExecutorAtResource):

    def __init__(self, job_id):
        self.host = None
        self.username = None
        self.known_hosts_filepath = None
        self.remote_path = os.path.join(ROOT_DIR, job_id)
        self.remote_client = None
        self.loop = asyncio.get_event_loop()

    # RESOURCE
    def set_resource(self, resource_params):  # coger y procesar el json
        self.host = resource_params["jm_location"]['host']
        self.username = resource_params["jm_credentials"]['username']
        self.known_hosts_filepath = resource_params["jm_credentials"]['known_hosts_filepath']



    def check(self):
        not_accessible = os.system(f"nc -z {self.host} 22")
        if not_accessible:  # if not_accessible is 0 the server is accessible else it is not
            return False
        else:
            return True

    def connect(self):
        self.remote_client = RemoteSSHClient(self.host, self.username, self.known_hosts_filepath, self.remote_path)
        self.loop.run_until_complete(self.remote_client.connect())
        return self.remote_client

    def disconnect(self):
        self.remote_client.disconnect()
        self.loop.close()

    # JOB EXECUTION
    def set_credentials(self, credentials):  # set username and password OR known_hosts
        """ Different from connecting to the resource, job submission may require identifying user, to check
            priority and the like
            """
        # TODO
        pass

    def get_quotas_for_current_credentials(self):
        # TODO
        pass

    def create_job_workspace(self, name):
        # the name is the job_id
        self.loop.run_until_complete(self.remote_client.make_directory(os.path.join(ROOT_DIR, name)))

    def remove_job_workspace(self, name):
        self.loop.run_until_complete(self.remote_client.remove_directory(os.path.join(ROOT_DIR, name)))

    def upload_file(self, job_context):
        i = job_context["transfer_state"]["idx"]
        local_path = job_context["process"]["inputs"]["data"][i]["path"]
        remote_path = job_context["process"]["inputs"]["data"][i]["remote_path"]
        return self.remote_client.upload_file(local_path, remote_path)

    def upload_directory(self, workspace, local_filename, remote_location):
        return self.remote_client.upload_directory(local_filename, remote_location)

    def move_file(self, remote_source, remote_destination):
        self.loop.run_until_complete(self.remote_client.move_file(remote_source, remote_destination))

    def remove_file(self, remote_filename):
        self.loop.run_until_complete(self.remote_client.remove_file(remote_filename))

    def submit(self, process):
        params = process["inputs"]["parameters"]
        return self.loop.run_until_complete(self.remote_client.run_client(params["script_file"], params["script_params"]))

    def job_status(self, native_id):
        return self.loop.run_until_complete(self.remote_client.command_status(native_id))

    def cancel_job(self, native_id):
        self.remote_client.kill_process(native_id)

    # SSH
    def download_file(self, job_context):
        i = job_context["transfer_state"]["idx"]
        local_path = job_context["process"]["outputs"][i]["path"]
        remote_path = job_context["process"]["outputs"][i]["remote_path"]
        self.loop.run_until_complete(self.remote_client.download_file(remote_path, local_path))

    def retrieve_directory(self, remote_dir, local_dir):
        self.loop.run_until_complete(self.remote_client.download_directory(remote_dir, local_dir))

    def exists(self, job_context):
        i = job_context["transfer_state"]["idx"]
        local_path = job_context["process"]["inputs"]["data"][i]["path"]
        remote_path = job_context["process"]["inputs"]["data"][i]["remote_path"]
        if os.path.exists(local_path):
            check = self.loop.run_until_complete(self.remote_client.exists_remotely(remote_path))

            if check:
                #the same_size method works ONLY with files and empty directories
                check &= self.loop.run_until_complete(
                    self.remote_client.same_size(local_path, remote_path))
            return check
        else:
            print(f"File {local_path} not found in your local system")
            return None

    def get_upload_files_list(self, job_context):
        return job_context["process"]["inputs"]["data"] +\
               job_context["process_params"]["parameters"]["script_files"]

    def get_download_files_list(self, job_context):
        #TODO
        pass

    def check_resource(self):
        """
        Can connect
        CPU and GPU use
        Available storage space

        :return:
        """
        # TODO:
        pass

    def prepare_execution_directory(self):
        # TODO: same as create_job_workspace?
        pass

    def check_executing_process(self):

        pass

    def remove_execution_directory(self, name):
        """
        Include parameters to keep intermediate files at some common directory (intermediate/<originating_job_id>/<file>
        :return:
        """
        # TODO: same as remove_job_workspace?
        pass