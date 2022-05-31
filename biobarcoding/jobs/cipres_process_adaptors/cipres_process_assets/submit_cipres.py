from sys import argv
import os
from python_cipres import client as cipres


def submit_cipres(base_url: str, username: str, password: str, app_id: str, app_name: str, vParams: {}, inputParams: {}):
    os.environ["URL"] = base_url
    os.environ["PASSWORD"] = password
    os.environ["KEY"] = app_id
    metadata = {"clientJobId": os.path.basename(os.path.dirname(__file__))}
    print("------------------ SUBMIT CIPRES ----------------------")
    print(inputParams)
    print(vParams)
    print(metadata)
    try:
        cipres_client = cipres.Client(app_name, app_id, username, password, base_url)
        cipres_client.verbose = True
        print("Cipres Client success")
    except Exception as e:
        print("Cipres Client")
        print(e)
        raise Exception()
    try:
        job_status = cipres_client.submitJob(vParams, inputParams, metadata)  # validateOnly=True
        job_status.show()
    except Exception as e:
        print("Submit Job")
        print(e)
    print("Writing to job_handle")
    with open("job_handle.txt", "w") as f:
        print(job_status.jobHandle)
        f.write(job_status.jobHandle)
    print("Succesful write")


if __name__ == "__main__":
    base_url, username, password, app_id, app_name, vParams, inputParams = argv[1:]
    submit_cipres(base_url, username, password, app_id, app_name, vParams, inputParams)
