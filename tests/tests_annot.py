import json

OBJECT_UUIDS = []


def identified(func):
    def wrapper(testful):
        login(testful)
        post_dataset(testful)
        func(testful)
        delete_dataset(testful)
        logout(testful)

    def login(testful):
        response = testful.put("/api/authn?user=test_user")
        while response.status_code >= 300:
            response = testful.put("/api/authn?user=test_user")

    def logout(testful):
        response = testful.delete("/api/authn")

    def post_dataset(testful):
        response = testful.post("/api/bos/sequences", json={'uniquename': 'test_seq_1', 'stock': 'test_seq_1'})
        response = testful.post("/api/bos/sequences", json={'uniquename': 'test_seq_2', 'stock': 'test_seq_2'})
        response = testful.get("/api/bos/sequences", json={'filter': {'uniquename': ('test_seq_1', 'test_seq_2')}})
        response = json.loads(response.data).get('content', {})
        global OBJECT_UUIDS
        OBJECT_UUIDS = [i.get('uuid') for i in response if i.get('uuid')]

    def delete_dataset(testful):
        response = testful.delete("/api/bos/sequences", json={'filter': {'uniquename': ('test_seq_1', 'test_seq_2')}})
        response = testful.delete("/api/individuals", json={'filter': {'uniquename': ('test_seq_1', 'test_seq_2')}})
        OBJECT_UUIDS.clear()

    return wrapper


def basic_cycle(api):
    assert api.get().status_code < 300
    assert api.get_one().status_code < 300
    assert api.create().status_code < 300
    assert api.get_one().status_code < 300
    assert api.modify().status_code < 300
    assert api.get_one().status_code < 300


def closing_cycle(api):
    assert api.remove().status_code < 300
    assert api.get_one().status_code < 300
    assert api.get().status_code < 300


def full_cycle(api):
    basic_cycle(api)
    closing_cycle(api)


@identified
def test_form_templates(testful):
    full_cycle(TemplatesAPI(testful))


@identified
def test_form_fields(testful):
    basic_cycle(TemplatesAPI(testful))
    full_cycle(FieldsAPI(testful))
    closing_cycle(TemplatesAPI(testful))


@identified
def test_annotations(testful):
    basic_cycle(TemplatesAPI(testful))
    basic_cycle(FieldsAPI(testful))
    full_cycle(AnnotationsAPI(testful))
    closing_cycle(FieldsAPI(testful))
    closing_cycle(TemplatesAPI(testful))


@identified
def test_form_relationships(testful):
    basic_cycle(TemplatesAPI(testful))
    basic_cycle(FieldsAPI(testful))
    full_cycle(FormRelationshipsAPI(testful))
    closing_cycle(FieldsAPI(testful))
    closing_cycle(TemplatesAPI(testful))


@identified
def test_relationships(testful):
    basic_cycle(TemplatesAPI(testful))
    basic_cycle(FieldsAPI(testful))
    basic_cycle(FormRelationshipsAPI(testful))
    full_cycle(RelationshipsAPI(testful))
    closing_cycle(FormRelationshipsAPI(testful))
    closing_cycle(FieldsAPI(testful))
    closing_cycle(TemplatesAPI(testful))


@identified
def test_end(testful):
    print('END')


class TemplatesAPI:

    def __init__(self, testful):
        self.testful = testful

    url = "/api/annotation_form_templates/"

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
        response = self.testful.post(self.url, data=self.post_req)
        print(response.data)
        return response

    def get_one(self):
        response = self.testful.get(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.data)
        return response

    def get(self):
        response = self.testful.get(self.url)
        return response

    def modify(self):
        response = self.testful.put(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}", data=self.put_req)
        print(response.data)
        return response

    def remove(self):
        response = self.testful.delete(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.data)
        return response


class FieldsAPI:

    def __init__(self, testful):
        self.testful = testful

    url = "/api/annotation_form_fields/"

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
        response = self.testful.post(self.url, data=self.post_req)
        print(response.data)
        return response

    def get_one(self):
        response = self.testful.get(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.data)
        return response

    def get(self):
        response = self.testful.get(self.url)
        return response

    def modify(self):
        response = self.testful.put(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}", data=self.put_req)
        print(response.data)
        return response

    def remove(self):
        response = self.testful.delete(f"{self.url}{self.post_req.get('db')}:{self.post_req.get('dbxref')}")
        print(response.data)
        return response


class AnnotationsAPI:

    def __init__(self, testful):
        self.testful = testful

    url = "/api/annotations/"

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
        "type": "text",
        "value": "maturase k",
    }

    put_req = {"value": "modified"}

    put_req_batch = [
        {
            "value": "<p>Hola, ¿qué tal?</p>",
            "name": "modified",
        },
        {
            "value": "<p>Hey</p>",
            "type": "text",
            "name": "append",
        },
    ]

    def create(self):
        # for payload in (self.template_post_req, self.field_post_req, self.post_req_text):
        #     response = self.testful.post(self.url, data=payload)
        #     print(response.data)
        data = [self.template_post_req, self.field_post_req, self.post_req_text]
        response = self.testful.post(self.url + OBJECT_UUIDS[0], data=str(data))
        print(response.data)
        return response

    def get_one(self):
        response = self.testful.get(f"{self.url}{OBJECT_UUIDS[0]}")
        print(response.data)
        return response

    def get(self):
        response = self.testful.get(self.url)
        return response

    def modify(self):

        response = self.testful.get(self.url)
        content = json.loads(response.data)['content']
        _id = content[-1].get('id') if content else ''

        response = self.testful.put(f"{self.url}{_id}",
                                    data=json.dumps({**self.put_req, **{'object_uuid': OBJECT_UUIDS}}))
        print('PUT: single modification')
        print(response.data)
        if response.status_code >= 300:
            return response

        response = self.testful.get(f"{self.url}{OBJECT_UUIDS[1]}")
        content = json.loads(response.data)['content']
        _id = content[-1].get('id') if content else ''
        response = self.testful.put(f"{self.url}{OBJECT_UUIDS[1]}",
                                    data=json.dumps([{**self.put_req_batch[0], **{'id': _id}}] + self.put_req_batch[1:]))
        print('PUT: modified and append')
        print(response.data)
        if response.status_code >= 300:
            return response

        return response

    def remove(self):
        response = self.testful.delete(f"{self.url}{OBJECT_UUIDS[0]}")
        print(response.data)
        response = self.testful.delete(f"{self.url}{OBJECT_UUIDS[1]}")
        print(response.data)
        return response


##
# Relationship annotations
##

class FormRelationshipsAPI:

    def __init__(self, testful):
        self.testful = testful

    url = "/api/annotation_form_relationships/"

    post_req = {
        "template": TemplatesAPI.post_req.get('name'),
        "field": FieldsAPI.post_req.get('name'),
        "name": "gene",
    }

    put_req = {
        "rev_name": "locus",
    }

    def create(self):
        response = self.testful.post(self.url, data=self.post_req)
        print(response.data)
        return response

    def get_one(self):
        response = self.testful.get(self.url)
        content = json.loads(response.data)['content']
        _id = content[-1].get('id') if content else ''

        response = self.testful.get(f"{self.url}{_id}")
        print(response.data)
        return response

    def get(self):
        response = self.testful.get(self.url)
        return response

    def modify(self):
        response = self.testful.get(self.url)
        content = json.loads(response.data)['content']
        _id = content[-1].get('id') if content else ''

        response = self.testful.put(f"{self.url}{_id}", data=self.put_req)
        print(response.data)
        return response

    def remove(self):
        response = self.testful.get(self.url)
        content = json.loads(response.data)['content']
        _id = content[-1].get('id') if content else ''

        response = self.testful.delete(f"{self.url}{_id}")
        print(response.data)
        return response


class RelationshipsAPI:

    def __init__(self, testful):
        self.testful = testful

    url = "/api/relationships/"

    put_req = json.dumps({
        "value": "modified",
    })

    def create(self):
        response = FormRelationshipsAPI(self.testful).get()
        content = json.loads(response.data)['content']
        _id = content[-1].get('id') if content else ''

        response = self.testful.post(self.url, data={"subject_uuid": OBJECT_UUIDS[0],
                                                     "object_uuid": OBJECT_UUIDS[1],
                                                     "form_rl_id": _id})
        print(response.data)
        return response

    def get_one(self):
        response = self.testful.get(self.url)
        content = json.loads(response.data)['content']
        _id = content[-1].get('id') if content else ''

        response = self.testful.get(f"{self.url}{_id}")
        print(response.data)
        return response

    def get(self):
        response = self.testful.get(self.url)
        return response

    def modify(self):
        response = self.testful.get(self.url)
        content = json.loads(response.data)['content']
        _id = content[-1].get('id') if content else ''

        response = self.testful.put(f"{self.url}{_id}", data=self.put_req)
        print(response.data)
        return response

    def remove(self):
        response = self.testful.get(self.url)
        content = json.loads(response.data)['content']
        _id = content[-1].get('id') if content else ''

        response = self.testful.delete(f"{self.url}{_id}")
        print(response.data)
        return response
