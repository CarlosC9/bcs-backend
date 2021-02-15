import psutil

from biobarcoding.jobs import JobExecutorAtResource
import os
import asyncio
import asyncssh


class RemoteSSHClient:
    """Client to interact with a remote host via SSH & SCP."""

    def __init__(self, host, username, known_hosts_filepath, remote_path):
        super()
        self.host = host
        self.username = username
        self.known_hosts_filepath = known_hosts_filepath
        self.remote_path = remote_path
        self.pid = None
        self.client = None
        self.sftp = None
        self.conn = None

    async def connect(self):
        """Open connection to remote host."""
        if self.conn is None:
            try:
                kwargs = {  # necesito known_hosts y los argumentos del constructor de RemoteSSHClient
                    'known_hosts': self.known_hosts_filepath,
                }
                self.conn, self.client = await asyncssh.create_connection(asyncssh.SSHClient, self.host, **kwargs)
                self.sftp = await self.conn.start_sftp_client()
                await self.sftp.chdir(self.remote_path)  # change remote working directory
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

    async def run_client(self, script_file):
        popen_pipe = os.popen(f"ssh {self.username}@{self.host} 'cd {self.remote_path} && chmod +x {script_file} && nohup ./{script_file} >/tmp/mtest2 </dev/null 2>/tmp/mtest2.err & echo $!'")
        pid = popen_pipe.readline().rstrip()
        print(f"PID: {pid}")
        return pid

    def command_status(self):
        """
        Get Job status
        :return status: -1: command not executed; 0: running; 1: not running
        """
        status = -1
        if self.pid is not None:
            status = os.system(f"ssh {self.username}@{self.host} 'ps --pid {self.pid}' >/dev/null 2>&1")
        else:
            print("No process being executing with the SSH Client")
        return status

    def kill_process(self):
        if self.pid is not None:
            os.system(f"ssh {self.username}@{self.host} 'kill -9 {self.pid}'")
            self.pid = None
        else:
            print("No command has been executed")

    '''async def upload_file(self, file_path, remote_file_path):
        """
        Upload multiple files to remote_path.

        :param file_path: List of paths to local files.
        :param remote_file_path: List of paths to local files.
        """
        if self.sftp is not None:
            await self.sftp.put(file_path, remotepath=os.path.join(self.remote_path, remote_file_path))
            print(f'Finished uploading {file_path} to {self.remote_path} on {self.host}')
        else:
            print("SSH connection not created")'''

    def upload_file_or_directory(self, local_path, remote_path):
            """
            Upload multiple files to remote_path.

            :param path: path to local folder or file. Folders must end with /
            :param remote_file_path: path to remote folder or file. Folders must end with /
            """
            cmd = f"nohup bash -c \"rsync {local_path} {self.username}@{self.host}:{os.path.join(self.remote_path, remote_path)}\" >/tmp/mtest2 </dev/null 2>/tmp/mtest2.err & echo $!"
            print(cmd)
            popen_pipe = os.popen(cmd)
            pid = popen_pipe.readline().rstrip()
            print(f"PID: {pid}")
            return pid

    async def upload_files_or_directories(self, files_dir_paths):
        """
        Upload multiple files and directories to remote_path.

        :param files_dir_paths: List of paths to local files or directories.
        :type files_dir_paths: List[str]
        """
        if self.sftp is not None:
            await self.sftp.put(files_dir_paths, recurse=True)
            print(
                f'Finished uploading {len(files_dir_paths)} files or directories to {self.remote_path} on {self.host}')
        else:
            print("SSH connection not created")

    async def download_file(self, remote_file, local_file):
        """Download file from remote host."""
        if self.sftp is not None:
            await self.sftp.get(os.path.join(self.remote_path, remote_file), localpath=local_file)
        else:
            print("SSH connection not created")

    async def download_files_or_directories(self, remote_file_dir, local_file_dir):
        """Download file from remote host."""
        if self.sftp is not None:
            await self.sftp.get(os.path.join(self.remote_path, remote_file_dir), localpath=local_file_dir, recurse=True)
        else:
            print("SSH connection not created")

    async def move_file(self, filepath, new_filepath):
        if self.sftp is not None:
            await self.sftp.rename(os.path.join(self.remote_path, filepath),
                                   os.path.join(self.remote_path, new_filepath))
        else:
            print("SSH connection not created")

    def remove_dir(self, dirname):
        popen_pipe = os.popen(
            f"nohup ssh {self.username}@{self.host} rm -r {dirname} >/tmp/mtest2 </dev/null 2>/tmp/mtest2.err & echo $!")
        pid = popen_pipe.readline().rstrip()
        print(f"PID: {pid}")
        return pid

    async def remove_file(self, file):
        if self.sftp is not None:
            await self.sftp.remove(os.path.join(self.remote_path, file))
        else:
            print("SSH connection not created")

    async def remove_files(self, files):
        """
        Upload multiple files to a remote directory.

        :param files: List of paths to local files.
        :type files: List[str]
        """
        if self.sftp is not None:
            statements = []
            for file in files:
                statements.append(self.remove_file(file))
            await asyncio.gather(*statements)
            print(f'Finished removing {len(files)} files from {self.remote_path} on {self.host}')
        else:
            print("SSH connection not created")

    async def check_file(self, file):
        if self.sftp is not None:
            return await self.sftp.isfile(os.path.join(self.remote_path, file))
        else:
            print("SSH connection not created")

    async def make_directory(self, dir_name):
        if self.sftp is not None:
            if not await self.exists_remotely(dir_name) and self.sftp.isdir(os.path.join(self.remote_path, dir_name)):
                # TODO: Se le puede asignar permisos, por ejemplo al working directory
                await self.sftp.mkdir(os.path.join(self.remote_path, dir_name))
            else:
                print("folder already exists")
        else:
            print("SSH connection not created")

    '''async def remove_directory(self, dir_path):
        if self.sftp is not None:
            await self.sftp.rmtree(dir_path)
        else:
            print("SSH connection not created")'''

    async def exists_remotely(self, name):
        return await self.sftp.exists(os.path.join(self.remote_path, name))

    async def same_size(self, local_name, remote_name):
        return os.path.getsize(local_name) == await self.sftp.getsize(os.path.join(self.remote_path, remote_name))


class JobExecutorWithSSH(JobExecutorAtResource):

    def __init__(self):
        self.host = None
        self.username = None
        self.known_hosts_filepath = None
        self.remote_path = None
        self.local_workspace = None
        self.remote_client = None
        self.up_files = None
        self.loop = asyncio.get_event_loop()

    # RESOURCE
    def set_resource(self, params):  # coger y procesar el json
        self.host = params['host']
        self.username = params['username']
        self.known_hosts_filepath = params['known_hosts_filepath']
        self.remote_path = params['remote_workspace']
        self.up_files = params['upload_files']
        self.local_workspace = params['local_workspace'] if 'local_workspace' in params else None

    def check(self):
        exit_status = os.system(f"nc -z {self.host} 22")
        if exit_status:  # if exit_status is 0 the server is accessible else it is not
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
        if self.remote_path == name:
            # it is empty because the remote_path is already a parameter in remote_client
            self.loop.run_until_complete(self.remote_client.make_directory(""))
            self.remote_client.upload_file_or_directory(self.local_workspace)
        else:
            print("The Job Workspace should be equal to remote_path")

    def remove_job_workspace(self, name):  # After Job is completed (or if Job was not started)
        if name != self.remote_path:
            raise Exception("JobExecutionWithSSH - upload_file - remote_location must be equal to self.remote_path")
            pid = self.remote_client.remove_dir(self.remote_path)
        return pid

    def upload_file(self, workspace, local_fileidx, remote_location):
        '''if remote_location != self.remote_path:
            raise Exception("JobExecutionWithSSH - upload_file - remote_location must be equal to self.remote_path")'''

        pid = self.remote_client.upload_file_or_directory(
                    self.up_files[local_fileidx]["local_path"], self.up_files[local_fileidx]["remote_path"])

        return pid

    def move_file(self, remote_source, remote_destination):
        self.loop.run_until_complete(self.remote_client.move_file(remote_source, remote_destination))

    def remove_file(self, remote_filename):
        self.loop.run_until_complete(self.remote_client.remove_file(remote_filename))

    def submit(self, workspace, params):
        return self.loop.run_until_complete(self.remote_client.run_client(params["script_file"]))

    '''def remote_job_status(self, native_id):
        return self.remote_client.command_status()'''

    def job_status(self, native_id):
        for proc in psutil.process_iter():
            try:
                if native_id == proc.pid:
                    return True
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                pass
        return False;

    def cancel_job(self, native_id):
        self.remote_client.kill_process()

    # SSH
    def send_files_or_directories(self, files):
        self.loop.run_until_complete(self.remote_client.upload_files_or_directories(self, files))

    def retrieve_files_or_directories(self, remote_files_dirs, local_dir):
        self.loop.run_until_complete(self.remote_client.download_files_or_directories(remote_files_dirs, local_dir))

    def remove_files(self, remote_files):
        self.loop.run_until_complete(self.remote_client.remove_files(remote_files))

    def check_files(self, remote_files):
        return self.loop.run_until_complete(self.remote_client.check_files(remote_files))

    def make_directories(self, dir_names):
        statements = []
        for dir_name in dir_names:
            statements.append(self.remote_client.make_directory(dir_name))
        self.loop.run_until_complete(asyncio.gather(*statements))

    def exists(self, file_idx):
        check = self.loop.run_until_complete(self.remote_client.exists_remotely(self.up_files[file_idx]["remote_path"]))
        if check:
            check &= self.loop.run_until_complete(self.remote_client.same_size(
                        self.up_files[file_idx]["local_path"], self.up_files[file_idx]["remote_path"]))
        return check

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

    def execute(self, commands):
        # TODO: same as submit?
        pass

    def check_executing_process(self):
        # TODO: same as create_job_workspace?
        pass

    def remove_execution_directory(self, name):
        """
        Include parameters to keep intermediate files at some common directory (intermediate/<originating_job_id>/<file>
        :return:
        """
        # TODO: same as remove_job_workspace?
        pass