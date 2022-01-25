from sys import argv
import os
import json
from python_cipres import client as cipres


def submit_cipres(base_url, username, password, app_id, app_name, vParams, inputParams):
    os.environ["URL"] = base_url
    os.environ["PASSWORD"] = password
    os.environ["KEY"] = app_id
    metadata = {}
    cipres_client = cipres.Client(app_name, app_id, username, password, base_url)
    job_status = cipres_client.submitJob(json.load(vParams), json.load(inputParams), metadata, validateOnly=True)  # TODO
    with open("job_handle.txt", "w") as f:
        f.write(job_status.jobHandle)


if __name__ == "__main__":
    base_url, username, password, app_id, app_name, vParams, inputParams = argv[1:]
    submit_cipres(base_url, username, password, app_id, app_name, vParams, inputParams)
