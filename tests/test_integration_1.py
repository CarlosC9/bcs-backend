import pytest

from biobarcoding.rest.main import create_app


def test_integration_1(testful):
    """
    Simulation of the sequence of calls GUI-Backend from load of GUI structure, data,
    to the launch of a dummy job and its results
    :return:
    """
    # TODO login()
    # TODO import_ontologies()
    # TODO import_taxonomy()
    # TODO import_sequences()
    # TODO get_menu()
    # TODO get_processes()
    # TODO get_process_input_form()
    # TODO get_resources_for_process()
    # submit_job(p, r, inputs)
    response = testful.post("/api/jobs/", json=dict(process="dummy", resource="local", inputs=dict()))
    # TODO wait until job completion
    # TODO read results
    print(response)


