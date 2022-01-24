import unittest
import requests

session = requests.Session()


def test_log_in():
    url = "http://localhost:5000/api/authn?user=test_user"

    response = session.put(url)
    print(response.text)

    return response


class TemplatesAPITest(unittest.TestCase):

    url = "http://localhost:5000/api/annotation_form_templates/"

    post_req = {
        # "attributes": {"tags": ["test", "testing"]},
        "name": "new_template",
        "cvterm_id": "",
        "dbxref_id": "",
        "db": "SO",
        "dbxref": "0000110",
        "object_type": '["sequence", "multiple-sequence-alignment"]',
    }

    put_req = {
        "description": "modified",
        "object_type": "sequence",
    }

    def test_create_template(self):
        response = session.post(self.url, data=self.post_req)
        print(response.text)

        self.assertEqual(response.status_code, 201)
        return None

    def test_get_template_by_id(self):
        response = session.get(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def test_get_templates(self):
        response = session.get(self.url)
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    def test_modify_template(self):
        response = session.put(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}", data=self.put_req)
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    # def test_detele_template(self):
    #     response = session.delete(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
    #     print(response.text)
    #
    #     self.assertEqual(response.status_code, 200)
    #     return None


class FieldsAPITest(unittest.TestCase):

    url = "http://localhost:5000/api/annotation_form_fields/"

    post_req = {
        "name": "new_field",
        "cvterm_id": "",
        "dbxref_id": "",
        "db": "SO",
        "dbxref": "0000001",
        "object_type": '["sequence", "multiple-sequence-alignment"]',
        "type": "tag",
        "range": "",
        "view_type": "select",
        "multiple": "0",
        "template": TemplatesAPITest.post_req.get('name'),
    }

    put_req = {
        "template": 'null',
        "description": "modified",
        "object_type": "sequence",
    }

    def test_create_field(self):
        response = session.post(self.url, data=self.post_req)
        print(response.text)

        self.assertEqual(response.status_code, 201)
        return None

    def test_get_field_by_id(self):
        response = session.get(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    def test_get_fields(self):
        response = session.get(self.url)
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    def test_modify_field(self):
        response = session.put(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}", data=self.put_req)
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    # def test_detele_field(self):
    #     response = session.delete(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
    #     print(response.text)
    #
    #     self.assertEqual(response.status_code, 200)
    #     return None


class AnnotationsAPITest(unittest.TestCase):

    url = "http://localhost:5000/api/annotations/"
    object_uuids = ["0bae8453-fdd5-42fa-9f89-6a1af2771183", "4a2b5bb0-452c-452e-b5c3-a00b7a644351"]

    post_req_template = {
        # "attributes": {"tags": ["test", "testing"]},
        # "object_uuid": object_uuids[0],
        "name": "new_annotation_template",
        "type": "template",
        "template": TemplatesAPITest.post_req.get('name'),
        "value": "matk",
    }

    post_req_field = {
        # "attributes": {"tags": ["test", "testing"]},
        # "object_uuid": object_uuids[0],
        "name": "new_annotation_field",
        "type": "field",
        "field": FieldsAPITest.post_req.get('name'),
        "value": "matk",
    }

    post_req_text = {
        # "attributes": {"tags": ["test", "testing"]},
        # "object_uuid": object_uuids[0],
        "name": "new_annotation_text",
        "type": "text",
        "value": "maturase k",
    }

    put_req = {
        "value": "modified",
        "object_uuid": str(object_uuids),
    }

    def test_create_annotations_by_uuid(self):

        data = [self.post_req_template, self.post_req_field, self.post_req_text]
        response = session.post(self.url + self.object_uuids[0], data=str(data))
        print(response.text)

        self.assertEqual(response.status_code, 201)
        return None

    # def test_create_annotation(self):
    #
    #     for payload in (self.post_req_template, self.post_req_field, self.post_req_text):
    #         response = session.post(self.url, data=payload)
    #         print(response.text)
    #
    #         self.assertEqual(response.status_code, 201)
    #
    #     return None

    def test_get_annotation_by_id(self):
        response = session.get(f"{self.url}{self.post_req_text.get('object_uuid')}")
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    def test_get_annotations(self):
        response = session.get(self.url)
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    def test_modify_annotation(self):
        response = session.put(f"{self.url}{self.post_req_text.get('object_uuid')}", data=self.put_req)
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    # def test_detele_annotation(self):
    #     response = session.delete(f"{self.url}{self.post_req_text.get('object_uuid')}")
    #     print(response.text)
    #
    #     self.assertEqual(response.status_code, 200)
    #     return None

