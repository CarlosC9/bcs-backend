import json
import unittest
import requests

session = requests.Session()
response = session.put("http://localhost:5000/api/authn?user=test_user")
while response.status_code >= 300:
    response = session.put("http://localhost:5000/api/authn?user=test_user")

OBJECT_UUIDS = ["652882bd-1cdc-4b94-a5ad-361fadedee28", "a3885d83-365c-41b8-a413-f6aeef14445d"]


class AnnotationSystemTest(unittest.TestCase):

    def basic_cycle(self, api):
        self.assertTrue(api.get().status_code < 300)
        self.assertTrue(api.get_by_id().status_code < 300)
        self.assertTrue(api.create().status_code < 300)
        self.assertTrue(api.get_by_id().status_code < 300)
        self.assertTrue(api.modify().status_code < 300)
        self.assertTrue(api.get_by_id().status_code < 300)

    def full_cycle(self, api):
        self.basic_cycle(api)
        self.assertTrue(api.remove().status_code < 300)
        self.assertTrue(api.get_by_id().status_code < 300)
        self.assertTrue(api.get().status_code < 300)

    def test_form_templates(self):
        self.full_cycle(TemplatesAPI())

    def test_form_fields(self):
        self.full_cycle(FieldsAPI())

    def test_annotations(self):
        self.basic_cycle(TemplatesAPI())
        self.basic_cycle(FieldsAPI())
        self.full_cycle(AnnotationsAPI())

    def test_form_relationships(self):
        self.basic_cycle(TemplatesAPI())
        self.basic_cycle(FieldsAPI())
        self.full_cycle(FormRelationshipsAPI())

    def test_relationships(self):
        self.basic_cycle(TemplatesAPI())
        self.basic_cycle(FieldsAPI())
        self.basic_cycle(FormRelationshipsAPI())
        self.full_cycle(RelationshipsAPI())


class TemplatesAPI:

    url = "http://localhost:5000/api/annotation_form_templates/"

    post_req = {
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

    def create(self):
        response = session.post(self.url, data=self.post_req)
        print(response.text)
        return response

    def get_by_id(self):
        response = session.get(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.text)
        return response

    def get(self):
        response = session.get(self.url)
        return response

    def modify(self):
        response = session.put(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}", data=self.put_req)
        print(response.text)
        return response

    def remove(self):
        response = session.delete(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.text)
        return response


class FieldsAPI:

    url = "http://localhost:5000/api/annotation_form_fields/"

    post_req = {
        "name": "new_field",
        "db": "SO",
        "dbxref": "0000001",
        "object_type": '["sequence", "multiple-sequence-alignment"]',
        "view_type": "select",
        "range": '["matk", "rbcl", "its"]'
    }

    put_req = {
        "template": TemplatesAPI.post_req.get('name'),
        "description": "modified",
        "object_type": "sequence",
    }

    def create(self):
        response = session.post(self.url, data=self.post_req)
        print(response.text)
        return response

    def get_by_id(self):
        response = session.get(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.text)
        return response

    def get(self):
        response = session.get(self.url)
        return response

    def modify(self):
        response = session.put(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}", data=self.put_req)
        print(response.text)
        return response

    def remove(self):
        response = session.delete(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.text)
        return response


class AnnotationsAPI:

    url = "http://localhost:5000/api/annotations/"

    template_post_req = {
        "name": "new_annotation_template",
        "template": TemplatesAPI.post_req.get('name'),
        "value": "matk",
    }

    field_post_req = {
        "name": "new_annotation_field",
        "field": FieldsAPI.post_req.get('name'),
        "value": "matk",
    }

    post_req_text = {
        "name": "new_annotation_text",
        "value": "maturase k",
    }

    put_req = json.dumps({
        "value": "modified",
        "object_uuid": OBJECT_UUIDS,
    })

    put_req_batch_1 = json.dumps([template_post_req, field_post_req, post_req_text])
    put_req_batch_2 = json.dumps([
        {
            "value": "<p>Hola, ¿qué tal?</p>",
            "name": "modified",
            "id": 1
        },
        {
            "value": "<p>Hey</p>",
            "name": "append",
        },
    ])

    def create(self):
        # for payload in (self.template_post_req, self.field_post_req, self.post_req_text):
        #     response = session.post(self.url, data=payload)
        #     print(response.text)
        data = [self.template_post_req, self.field_post_req, self.post_req_text]
        response = session.post(self.url + OBJECT_UUIDS[0], data=str(data))
        print(response.text)
        return response

    def get_by_id(self):
        response = session.get(f"{self.url}{OBJECT_UUIDS[0]}")
        print(response.text)
        return response

    def get(self):
        response = session.get(self.url)
        return response

    def modify(self):

        response = session.put(f"{self.url}1", data=self.put_req)
        print('PUT: single modification')
        print(response.text)
        if response.status_code >= 300:
            return response

        response = session.put(f"{self.url}{OBJECT_UUIDS[1]}", data=self.put_req_batch_1)
        print('PUT: create')
        print(response.text)
        if response.status_code >= 300:
            return response

        response = session.put(f"{self.url}{OBJECT_UUIDS[1]}", data=self.put_req_batch_2)
        print('PUT: modified and append')
        print(response.text)
        if response.status_code >= 300:
            return response

        return response

    def remove(self):
        response = session.delete(f"{self.url}{OBJECT_UUIDS[0]}")
        print(response.text)
        return response


##
# Relationship annotations
##

class FormRelationshipsAPI:

    url = "http://localhost:5000/api/annotation_form_relationships/"

    post_req = {
        "template": "new_template",
        "field": "new_field",
        "name": "gene",
    }

    put_req = {
        "rev_name": "locus",
    }

    def create(self):
        response = session.post(self.url, data=self.post_req)
        print(response.text)
        return response

    def get_by_id(self):
        response = session.get(self.url)
        content = json.loads(response.text)['content']
        _id = content[-1].get('id')

        response = session.get(f"{self.url}{_id}")
        print(response.text)
        return response

    def get(self):
        response = session.get(self.url)
        return response

    def modify(self):
        response = session.get(self.url)
        content = json.loads(response.text)['content']
        _id = content[-1].get('id')

        response = session.put(f"{self.url}{_id}", data=self.put_req)
        print(response.text)
        return response

    def remove(self):
        response = session.get(self.url)
        content = json.loads(response.text)['content']
        _id = content[-1].get('id')

        response = session.delete(f"{self.url}{_id}")
        print(response.text)
        return response


class RelationshipsAPI:

    url = "http://localhost:5000/api/relationships/"

    post_req = {
        "name": "new_relationship",
        "template": TemplatesAPI.post_req.get('name'),
        "value": "matk",
    }

    put_req = json.dumps({
        "value": "modified",
        "object_uuid": OBJECT_UUIDS,
    })

    put_req_batch_1 = json.dumps([post_req])
    put_req_batch_2 = json.dumps([
        {
            "value": "<p>Hola, ¿qué tal?</p>",
            "name": "modified",
            "id": 1
        },
        {
            "value": "<p>Hey</p>",
            "name": "append",
        },
    ])

    def create_by_uuid(self):

        data = [self.post_req]
        response = session.post(self.url + OBJECT_UUIDS[0], data=str(data))
        print(response.text)
        return response

    def create(self):
        for payload in (self.post_req, FieldsAPI.post_req, self.post_req_text):
            response = session.post(self.url, data=payload)
            print(response.text)
        return response

    def get_by_id(self):
        response = session.get(f"{self.url}{OBJECT_UUIDS[0]}")
        print(response.text)
        return response

    def get(self):
        response = session.get(self.url)
        return response

    def modify(self):
        response = session.put(f"{self.url}1", data=self.put_req)
        print(response.text)
        if response.status_code >= 300:
            return response

        response = session.put(f"{self.url}{OBJECT_UUIDS[1]}", data=self.put_req_batch_1)
        print('PUT: create')
        print(response.text)
        if response.status_code >= 300:
            return response

        response = session.put(f"{self.url}{OBJECT_UUIDS[1]}", data=self.put_req_batch_2)
        print('PUT: modified and append')
        print(response.text)
        return response

    def remove(self):
        response = session.get(self.url)
        content = json.loads(response.text)['content']
        _id = content[-1].get('id')

        response = session.delete(f"{self.url}{_id}")
        print(response.text)
        return response