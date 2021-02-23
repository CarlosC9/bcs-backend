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
        self.job_status_dir = os.path.join("/tmp", username + "_jobs_status")
        self.client = None
        self.sftp = None
        self.conn = None
        self.last_job_remotely = None

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
                if not os.path.isdir(self.job_status_dir): #create folder for saving local jobs exit status
                    os.mkdir(self.job_status_dir)
                    #TODO meter exit status en remote path directamente
                if not await self.exists_remotely(self.job_status_dir): #create folder for saving remote jobs exit status
                    await self.make_directory(self.job_status_dir)
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
        #TODO cambiar ruta mtest2 y mtest2.err a self.remote_path y mirar como dar parámetros al script_file de manera standard
        cmd = f"ssh {self.username}@{self.host} 'cd {self.remote_path} && chmod +x {script_file} && (nohup ./{script_file}  >/tmp/mtest2 </dev/null 2>/tmp/mtest2.err & echo $!; wait $!; echo $? >> {self.job_status_dir}/$!.exit_status)'"
        print(cmd)
        popen_pipe = os.popen(cmd)
        pid = popen_pipe.readline().rstrip()
        self.last_job_remotely = True
        print(f"PID: {pid}")
        return pid

    async def command_status(self, pid):
        """
        Get Job status
        :return exit_status
        """
        exit_status = -1#means no process with this pid
        if pid is not None:
            if self.last_job_remotely:
                #TODO el directorio de job status va a ser el mismo self.remote_path
                if await self.exists_remotely(f"{self.job_status_dir}/{pid}.exit_status"):
                    cmd = f"ssh {self.username}@{self.host} 'cat {self.job_status_dir}/{pid}.exit_status'"
                    popen_pipe = os.popen(cmd)
                    exit_status = popen_pipe.readline().strip()
                else:
                    exit_status = -2 #-2 means running
            else:
                if os.path.isfile(f"{self.job_status_dir}/{pid}.exit_status"):
                    cmd = f"cat {self.job_status_dir}/{pid}.exit_status"
                    popen_pipe = os.popen(cmd)
                    exit_status = popen_pipe.readline().strip()
                else:
                    exit_status = -2#-2 means running
            print(f"Exit Status: {exit_status}")
        else:
            print(f"No process being executing with PID: {pid}")
        return int(exit_status)

    def kill_process(self, pid):
        if pid is not None:
            if self.last_job_remotely:
                os.system(f"ssh {self.username}@{self.host} 'kill -9 {pid}'")
            else:
                os.system(f"kill -9 {pid}")
        else:
            print("No command has been executed")

    def upload_directory(self, local_path, remote_path):
            """
            Upload multiple files to remote_path.

            :param path: path to local folder.
            :param remote_file_path: path to remote folder relative to self.remote_path.
            """
            remote_path = os.path.join(self.remote_path, remote_path, "")#ensure it always ends with / (the separator)
            remote_dir = os.path.join(self.remote_path, os.path.split(remote_path)[0])
            if remote_dir != "" and False:
                create_remote_dir_cmd = f"--rsync-path='mkdir -p {remote_dir} & rsync'"
            else:
                create_remote_dir_cmd = ""
            #-av for folders TODO mirar si es necesaria la v
            cmd = f"(nohup bash -c \"rsync -av {create_remote_dir_cmd} {local_path} {self.username}@{self.host}:{remote_path}\" >/tmp/mtest2 </dev/null 2>/tmp/mtest2.err & echo $!; wait $!; echo $? >> {self.job_status_dir}/$!.exit_status)"
            print(cmd)
            popen_pipe = os.popen(cmd)
            self.last_job_remotely = False
            pid = popen_pipe.readline().rstrip()
            print(f"PID: {pid}")
            return pid

    async def download_directory(self, remote_dir, local_dir):
        """Download file from remote host."""
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

    def upload_file(self, local_path, remote_path):
            """
            Upload file to remote_path.

            :param path: path to local file.
            :param remote_file_path: path to file relative to self.remote_path.
            """
            remote_path = os.path.join(self.remote_path, remote_path)
            remote_dir = os.path.join(self.remote_path, os.path.split(remote_path)[0])
            if remote_dir != "" and False:
                create_remote_dir_cmd = f"--rsync-path='mkdir -p {remote_dir} & rsync'"
            else:
                create_remote_dir_cmd = ""
            #-av for folders
            cmd = f"(nohup bash -c \"rsync {create_remote_dir_cmd} {local_path} {self.username}@{self.host}:{remote_path}\" >/tmp/mtest2 </dev/null 2>/tmp/mtest2.err & echo $!; wait $!; echo $? >> {self.job_status_dir}/$!.exit_status)"
            print(cmd)
            popen_pipe = os.popen(cmd)
            self.last_job_remotely = False
            pid = popen_pipe.readline().rstrip()
            print(f"PID: {pid}")
            return pid

    async def download_file(self, remote_file, local_file):
        """Download file from remote host."""
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
        if self.sftp is not None:
            await self.sftp.rename(os.path.join(self.remote_path, filepath),
                                   os.path.join(self.remote_path, new_filepath))
        else:
            print("SSH connection not created")

    async def remove_file(self, file):
        if self.sftp is not None:
            await self.sftp.remove(os.path.join(self.remote_path, file))
        else:
            print("SSH connection not created")

    async def check_file(self, file):
        if self.sftp is not None:
            return await self.sftp.isfile(os.path.join(self.remote_path, file))
        else:
            print("SSH connection not created")

    async def make_directory(self, dir_name):
        if self.sftp is not None:
            if not await self.sftp.isdir(os.path.join(self.remote_path, dir_name)):
                # TODO: Se le puede asignar permisos, por ejemplo al working directory
                await self.sftp.mkdir(os.path.join(self.remote_path, dir_name))
            else:
                print("folder already exists")
        else:
            print("SSH connection not created")

    async def remove_directory(self, dir_path):
        if self.sftp is not None:
            await self.sftp.rmtree(os.path.join(self.remote_path, dir_path))
        else:
            print("SSH connection not created")

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
        self.loop.run_until_complete(self.remote_client.make_directory(name))

    def remove_job_workspace(self, name):
        self.loop.run_until_complete(self.remote_client.remove_directory(name))

    def upload_file(self, workspace, local_filename, remote_location):
        return self.remote_client.upload_file(local_filename, remote_location)

    def upload_directory(self, workspace, local_filename, remote_location):
        return self.remote_client.upload_directory(local_filename, remote_location)

    def move_file(self, remote_source, remote_destination):
        self.loop.run_until_complete(self.remote_client.move_file(remote_source, remote_destination))

    def remove_file(self, remote_filename):
        self.loop.run_until_complete(self.remote_client.remove_file(remote_filename))

    def submit(self, workspace, params):
        return self.loop.run_until_complete(self.remote_client.run_client(params["script_file"]))

    def job_status(self, native_id):
        return self.loop.run_until_complete(self.remote_client.command_status(native_id))

    def cancel_job(self, native_id):
        self.remote_client.kill_process(native_id)

    # SSH
    def retrieve_file(self, remote_file, local_file):
        self.loop.run_until_complete(self.remote_client.download_file(remote_file, local_file))

    def retrieve_directory(self, remote_dir, local_dir):
        self.loop.run_until_complete(self.remote_client.download_directory(remote_dir, local_dir))

    def exists(self, local_path, remote_path):
        check = self.loop.run_until_complete(self.remote_client.exists_remotely(remote_path))
        #todo mirar que
        if check:
            check &= self.loop.run_until_complete(self.remote_client.same_size(local_path, remote_path))
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

    def check_executing_process(self):

        pass

    def remove_execution_directory(self, name):
        """
        Include parameters to keep intermediate files at some common directory (intermediate/<originating_job_id>/<file>
        :return:
        """
        # TODO: same as remove_job_workspace?
        pass