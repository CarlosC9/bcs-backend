from ssh_resource import RemoteClient
import os

remote_working_dir = "/home/dreyes/test_data"

remote_client = ("balder", "dreyes", "/home/daniel/.ssh/id_rsa", remote_working_dir)
dirname = os.path.join(os.path.dirname(__file__), "../../biobarcoding/jobs/test_data")

files = [os.path.join(dirname, f) for f in os.listdir(dirname) if os.path.isfile(os.path.join(dirname, f))]

commands = [f'python3 {remote_working_dir}/count_lines.py',
            f'python3 {remote_working_dir}/add_line.py',
            f'python3 {remote_working_dir}/count_lines.py']

remote_client.execute_commands(commands)
remote_client.download_file("myfile.txt", os.path.join(os.path.dirname(__file__), "../../biobarcoding/jobs/test_data", "newmyfile.txt"))
remote_client.remove_directory(remote_working_dir)
