import pytest

from biobarcoding.rest.main import create_app


def test_integration_1(testful):
    """
    Pseudocode of an scenario of the main use case,
    imitating the sequence of calls GUI->Backend
    from login, to the launch of a dummy job and its results

    Pseudocode of an scenario of the main use case
    :param testful:
    :return:
    """
    # "login" (returns an identity)
    response = testful.put("/api/authn?user=test_user")
    # Check if the response is as expected result
    assert response.status_code == 200
    assert response.json["status"] == 'success'
    # Call "am i logged" (pass headers; returns a boolean)
    response = testful.get("/api/authn")
    assert response.status_code == 200
    assert response.json["result"] == 'True'
    # Call "get ontologies" (returns a list)
    response = testful.get("/api/ontologies/")
    assert response.status_code == 200
    assert type(response.json) == list
    # NOTE: Call "import ontologies" not published as an API. It is directly performed with scripts, at deploy
    # Call "get sequences" (returns a list of codes, JSON format by default)
    # TODO How many should be returned (there must be a limit; look for pagination support)
    # TODO Also, prepare to filter on Identity: each sequence has to be checked for "list" permission
    response = testful.get("/api/bos/sequences/")
    assert response.status_code == 200
    assert type(response.json) == list
    # Detail of a sequence, JSON format
    sid = None
    response = testful.get(f"/api/bos/sequences/{sid}")
    assert response.status == 200
    assert type(response.json) == dict
    # Call "import sequences" (import a FASTA file)
    # TODO For the upload mechanism, the same used in Ontology
    # TODO Optionally import them with metadata like Collection, Identity (importer)
    # TODO Update already existing sequences. This is very important to
    response = testful.post("/api/bos/sequences/")
    assert response.status == 200
    # Call "get processes" (return a list of processes in the platform, filtered by Identity permissions)
    response = testful.get("/api/processes/")
    assert response.status_code == 200
    assert type(response.json) == list
    # Call "get process 'i' input form" (return the schema)
    # TODO The input form should be in JSONSchema format or similar
    pid = 10
    response = testful.get(f"/api/processes/{pid}")
    assert response.status_code == 200
    assert type(response.json) == dict
    # Call "get resources" (return resources supporting that process)
    response = testful.get(f"/api/processes/{pid}/resources/")
    assert response.status_code == 200
    assert type(response.json) == list
    # Call "submit asynchronous process" (returns Job ID)
    rid = None
    # TODO Because it is MSA, a list of seq IDs could be one of the parameters.
    #  Another possibility could be the filter to select the sequences
    # TODO     - Create Job Object and submit to Celery which will send to Galaxy
    #        - Export sequences as FASTA
    #        - Upload to Galaxy Server
    #        - Call process with gathered parameters
    #        - Administrative task in Celery to check currently running Jobs, to put in "ready" mode
    response = testful.post("/api/jobs", json=dict(process=pid, resource=rid, params={}))
    assert response.status == 200
    assert type(response.json) == dict
    # Call "get jobs". The default filter would be "not finished or cancelled jobs, visible by the user"
    response = testful.get("/api/jobs")
    assert response.status == 200
    assert type(response.json) == list
    # Call "get job ID"
    jid = None
    response = testful.get(f"/api/jobs/{jid}")
    assert response.status == 200
    assert type(response.json) == dict
    # Wait until finished
    while True:
        response = testful.get(f"/api/jobs/{jid}")
        job_status = response.json["status"]
        if job_status in ["succesfully_finished", "succesfully_finished", "cancelled"]:
            break
        import time
        time.sleep(2)
    # Call "read job log" (returns object IDs and their types; and the log)
    response = testful.get(f"/api/jobs/{jid}/results")
    assert response.status == 200
    assert type(response.json) == dict
    # Call "logout" (returns a boolean)
    response = testful.delete("/api/authn")
    # Check if the response is as expected result
    assert response.status == 200
    assert response.json["result"] == 'True'
