import unittest
import requests
import json


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
        "object_type": ["sequence", "multiple-sequence-alignment"],
    }

    put_req = {
        "description": "modified",
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

    def test_detele_template(self):
        return None
        response = session.delete(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None


class FieldsAPITest(unittest.TestCase):

    url = "http://localhost:5000/api/annotation_form_fields/"

    post_req = {
        "name": "new_field",
        "cvterm_id": "",
        "dbxref_id": "",
        "db": "SO",
        "dbxref": "0000001",
        "object_type": ["sequences", "multiple-sequence-alignments"],
        "type": "tag",
        "range": "tag",
        "view_type": "tag",
        "multiple": "0",
        "template_id": "1",
    }

    put_req = {
        "description": "modified",
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

    def test_detele_field(self):
        return None
        response = session.delete(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None


class AnnotationsAPITest(unittest.TestCase):

    url = "http://localhost:5000/api/annotations/"

    post_req_template = {
        # "attributes": {"tags": ["test", "testing"]},
        "object_uuid": "99856669-d9d9-45f6-beba-24e0a64a669e",
        "name": "new_annotation_template",
        "type": "template",
        "db": "SO",
        "dbxref": "0000110",
        "value": "matk",
    }

    post_req_field = {
        # "attributes": {"tags": ["test", "testing"]},
        "object_uuid": "99856669-d9d9-45f6-beba-24e0a64a669e",
        "name": "new_annotation_field",
        "type": "field",
        "db": "SO",
        "dbxref": "0000001",
        "value": "matk",
    }

    post_req_text = {
        # "attributes": {"tags": ["test", "testing"]},
        "object_uuid": "99856669-d9d9-45f6-beba-24e0a64a669e",
        "name": "new_annotation_text",
        "type": "text",
        "value": "maturase k",
    }

    put_req = {
        "rank": "modified",
    }

    def test_create_annotation(self):

        for payload in (self.post_req_template, self.post_req_field, self.post_req_text):
            response = session.post(self.url, data=payload)
            print(response.text)

            self.assertEqual(response.status_code, 201)

        return None

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

    def test_detele_annotation(self):
        return None
        response = session.delete(f"{self.url}{self.post_req_text.get('object_uuid')}")
        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None


if __name__ == '__main__':
    TemplatesAPITest.test_detele_template()
    TemplatesAPITest.test_create_template()
    TemplatesAPITest.test_get_templates()
    TemplatesAPITest.test_get_template_by_id()
    TemplatesAPITest.test_modify_template()
    TemplatesAPITest.test_detele_template()

    FieldsAPITest.test_detele_field()
    FieldsAPITest.test_create_field()
    FieldsAPITest.test_get_fields()
    FieldsAPITest.test_get_field_by_id()
    FieldsAPITest.test_modify_field()
    FieldsAPITest.test_detele_field()

    AnnotationsAPITest.test_detele_annotation()
    AnnotationsAPITest.test_create_annotation()
    AnnotationsAPITest.test_get_annotations()
    AnnotationsAPITest.test_get_annotation_by_id()
    AnnotationsAPITest.test_modify_annotation()
    AnnotationsAPITest.test_detele_annotation()
