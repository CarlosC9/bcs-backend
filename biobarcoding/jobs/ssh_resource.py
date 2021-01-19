<<<<<<< edf873c1026a63e6396c575774a2c944bddaed79
def send_files():
    pass
=======
from biobarcoding.jobs import JobExecutorAtResource
import os
from paramiko import SSHClient, AutoAddPolicy, RSAKey, SFTPClient
from paramiko.auth_handler import AuthenticationException, SSHException
from scp import SCPClient, SCPException
from stat import S_ISDIR

>>>>>>> ssh constructing


<<<<<<< edf873c1026a63e6396c575774a2c944bddaed79
def retrieve_files():
    pass
=======
    def __init__(self, host, user, ssh_key_filepath, remote_path):
        self.host = host
        self.user = user
        self.ssh_key_filepath = ssh_key_filepath
        self.remote_path = remote_path
        self.client = None
        self.scp = None
        self.conn = None
        self._upload_ssh_key()
    
    def _get_ssh_key(self):
        """
        Fetch locally stored SSH key.
        """
        try:
            self.ssh_key = RSAKey.from_private_key_file(self.ssh_key_filepath)
            #logger.info(f'Found SSH key at self {self.ssh_key_filepath}')
            print(f'Found SSH key at self {self.ssh_key_filepath}')
        except SSHException as error:
            #logger.error(error)
            print(error)
>>>>>>> ssh constructing


<<<<<<< edf873c1026a63e6396c575774a2c944bddaed79
def remove_files():
    pass

=======
    #TODO: BORRAR
    def _upload_ssh_key(self):
        try:
            os.system(f'ssh-copy-id -i {self.ssh_key_filepath} {self.user}@{self.host}>/dev/null 2>&1')
            os.system(f'ssh-copy-id -i {self.ssh_key_filepath}.pub {self.user}@{self.host}>/dev/null 2>&1')
            #logger.info(f'{self.ssh_key_filepath} uploaded to {self.host}')
            print(f'{self.ssh_key_filepath} uploaded to {self.host}')
        except FileNotFoundError as error:
            #logger.error(error)
            print(error)

    def connect(self):
        """Open connection to remote host."""
        if self.conn is None:
            try:
                self.client = SSHClient()
                self.client.load_system_host_keys()
                self.client.set_missing_host_key_policy(AutoAddPolicy())
                self.client.connect(
                    self.host,
                    username=self.user,
                    key_filename=self.ssh_key_filepath,
                    look_for_keys=True,
                    timeout=5000
                )
                self.scp = SCPClient(self.client.get_transport())
            except AuthenticationException as error:
                #logger.error(f'Authentication failed: \
                #    did you remember to create an SSH key? {error}')
                raise error
        return self.client
>>>>>>> ssh constructing

def check_files():
    pass


<<<<<<< edf873c1026a63e6396c575774a2c944bddaed79
def check_resource():
    """
    Can connect
    CPU and GPU use
    Available storage space

    :return:
    """
=======
        :param commands: List of unix commands as strings.
        :type commands: List[str]
        """
        #self.conn = self._connect()
        if self.client is not None:
            for cmd in commands:
                stdin, stdout, stderr = self.client.exec_command(cmd)
                stdout.channel.recv_exit_status()
                errors = stderr.readlines()
                for line in errors:
                    print(f'INPUT: {cmd} | ERROR: {line}')
                response = stdout.readlines()
                for line in response:
                    print(f'INPUT: {cmd} | OUTPUT: {line}')
        else:
            print("SSH connection not created")

    def _upload_single_file(self, file):
        """Upload a single file to a remote directory."""
        if self.scp is not None:
            upload = None
            try:
                self.scp.put(
                    file,
                    recursive=True,
                    remote_path=self.remote_path
                )
                upload = file
            except SCPException as error:
                print(error)
                raise error
            finally:
                print(f'Uploaded {file} to {self.remote_path}')
                return upload
        else:
            print("SSH connection not created")
>>>>>>> ssh constructing


<<<<<<< edf873c1026a63e6396c575774a2c944bddaed79
def prepare_execution_directory():
    pass


def execute():
    pass


def check_executing_process():
    pass


def remove_execution_directory():
    """
    Include parameters to keep intermediate files at some common directory (intermediate/<originating_job_id>/<file>
    :return:
    """
    pass

=======
        :param files: List of paths to local files.
        :type files: List[str]
        """
        #self.conn = self._connect()
        uploads = [self._upload_single_file(file) for file in files]
        #logger.info(f'Finished uploading {len(uploads)} files to {self.remote_path} on {self.host}')
        print(f'Finished uploading {len(uploads)} files to {self.remote_path} on {self.host}')

    def download_file(self, remote_file, local_file):
        """Download file from remote host."""
        if self.scp is not None:
            self.scp.get(os.path.join(self.remote_path, remote_file), local_path=local_file)
        else:
            print("SSH connection not created")

    #TODO: lo borro?
    def move_file(self, file, new_dir):
        if self.client is not None:
            sftp = self.client.open_sftp()
            sftp.rename(os.path.join(self.remote_path, file), os.path.join(new_dir, file))
            sftp.close()
        else:
            print("SSH connection not created")

    def remove_file(self, file):
        if self.client is not None:
            sftp = self.client.open_sftp()
            sftp.remove(os.path.join(self.remote_path, file))
            sftp.close()
        else:
            print("SSH connection not created")

    def bulk_remove(self, files):
        """
        Upload multiple files to a remote directory.

        :param files: List of paths to local files.
        :type files: List[str]
        """
        if self.client is not None:
            removed_files = [self.remove_file(file) for file in files]
            print(f'Finished removing {len(removed_files)} files from {self.remote_path} on {self.host}')
        else:
            print("SSH connection not created")

    def check_file(self, file):
        if self.client is not None:
            sftp = self.client.open_sftp()
            try:
                print(sftp.stat(os.path.join(self.remote_path, file)))
                return True
            except IOError:
                return False
            sftp.close()
        else:
            print("SSH connection not created")

    def make_directory(self, dir_name):
        if self.client is not None:
            sftp = self.client.open_sftp()
            sftp.mkdir(os.path.join(self.remote_path, dir_name))
            sftp.close()
        else:
            print("SSH connection not created")

    def _remove_directory(self, dir_name, sftp):
        files = sftp.listdir_attr(dir_name)

        for f in files:
            mode = f.st_mode
            if S_ISDIR(mode):
                self._remove_directory(os.path.join(dir_name, f.filename))
            else:
                sftp.remove(os.path.join(dir_name, f.filename))

        sftp.rmdir(dir_name)

    #TODO: probar bien
    def remove_directory(self, dir_name):
        if self.client is not None:
            sftp = self.client.open_sftp()
            path = os.path.join(self.remote_path, dir_name)
            self._remove_directory(path, sftp)
            sftp.close()
        else:
            print("SSH connection not created")
'''
class JobExecutionWithSSH(JobExecutorAtResource):

    #RESOURCE
    def set_resource(self, params):
        pass

    def check(self):
        pass
    #TODO: quitar parametros?
    def connect(self, host, user, ssh_key_filepath, remote_path):
        remote_client = RemoteClient(host, user, ssh_key_filepath, remote_path)
        remote_client.connect()
        return remote_client

    def disconnect(self, remote_client):
        remote_client.disconnect()

    # JOB EXECUTION
    def set_credentials(self, credentials):
        """ Different from connecting to the resource, job submission may require identifying user, to check
            priority and the like
            """
        pass

    def get_quotas_for_current_credentials(self):
        pass

    def create_job_workspace(self, name):
        remote_client = self.connect(#TODO: params)
        remote_client.make_directory(name)
        self.disconnect()

    def remove_job_workspace(self, name):  # After Job is completed (or if Job was not started)
        remote_client = self.connect(#TODO: params)
        remote_client.remove_directory(name)
        self.disconnect(remote_client)

    def upload_file(self, workspace, local_filename, remote_location):
        remote_client = self.connect(_, _, _, os.path.join(workspace, remote_location))
        remote_client.bulk_upload([local_filename])
        self.disconnect(remote_client)

    def move_file(self, remote_source, remote_destination):
        remote_client = self.connect(#TODO: params)
        remote_client.move_file(remote_source, remote_destination)
        self.disconnect(remote_client)

    def remove_file(self, remote_filename):
        remote_client = self.connect(#TODO: params)
        remote_client.remove_directory(remote_filename)
        self.disconnect(remote_client)


    def submit(self, workspace, params):
        #TODO: depende de como queramos los params

    def job_status(self, native_id):
        pass

    def cancel_job(self, native_id):
        remote_client = self.connect(_, _, user, _)
        #TODO: con docker usaremos el id
        exec_command("killall -u %s tail" % user)
        self.disconnect(remote_client)


    ## SSH

    def send_files(self, files):
        remote_client = self.connect(#TODO: params)
        remote_client.bulk_upload(self, files)
        self.disconnect(remote_client)


    def retrieve_files(self, remote_files):
        remote_client = self.connect(#TODO: params)
        for f in files:
            remote_client.download_file(self, f)
        self.disconnect(remote_client)


    def remove_files(self, remote_files):
        remote_client = self.connect(#TODO: params)
        for f in files:
            remote_client.remove_file(self, f)
        self.disconnect(remote_client)


    def check_files(self, remote_files):
        remote_client = self.connect(#TODO: params)
        for f in files:
            remote_client.remove_file(self, f)
        self.disconnect(remote_client)


    def check_resource():
        """
        Can connect
        CPU and GPU use
        Available storage space

        :return:
        """
        pass


    def prepare_execution_directory():
        pass


    def execute(self, commands):
        remote_client = self.connect(#TODO: params)
        remote_client.execute_commands(commands)
        self.disconnect(remote_client)


    def check_executing_process():
        pass


    def remove_execution_directory(self, name):
        """
        Include parameters to keep intermediate files at some common directory (intermediate/<originating_job_id>/<file>
        :return:
        """

        pass
'''
>>>>>>> ssh constructing
