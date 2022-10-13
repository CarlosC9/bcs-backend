import collections
import json
import os
from typing import List, Tuple, Union
import uuid

import sqlalchemy
from multidict import MultiDict, CIMultiDict
from sqlalchemy import and_
from sqlalchemy.pool import QueuePool

from .. import get_global_configuration_variable, engine, case_sensitive
from . import ROOT
from ..db_models import ORMBase
from ..db_models.jobs import StatusChecker
from ..db_models.jobs import ComputeResource, JobManagementType, Process, ProcessInComputeResource
# #####################################################################################################################
# >>>> DATABASE FUNCTIONS <<<<
# #####################################################################################################################


RESOURCE_PROCESSES_DICT = {
    "localhost - galaxy": [
    ],
    "balder - galaxy": [
    ],
    "localhost - ssh": [
        "MSA ClustalW",
        "PAUP Parsimony",
        "Phylogenetic Diversity Analyzer",
        "MSA ClustalW + PAUP Parsimony",
        "PAUP Parsimony + Phylogenetic Diversity Analyzer",
        "MSA ClustalW + PAUP Parsimony + Phylogenetic Diversity Analyzer"
    ],
    "balder - ssh": [
        "MSA ClustalW",
        "PAUP Parsimony",
        "Phylogenetic Diversity Analyzer",
        "MSA ClustalW + PAUP Parsimony",
        "PAUP Parsimony + Phylogenetic Diversity Analyzer",
        "MSA ClustalW + PAUP Parsimony + Phylogenetic Diversity Analyzer"
    ],
    "balder - slurm": [
        "MAFFT",
        "Mr Bayes",
        "BLAST",
        "Beast"
    ],
    "Teide HPC - slurm": [
        "MAFFT",
        "Mr Bayes",
        "BLAST",
        "Beast"
    ],
    "San Diego Supercomputer Center - cipres": [
        "Mr Bayes",
        "MAFFT"
    ]
}

PROCESSES_INPUTS = {
    "c8df0c20-9cd5-499b-92d4-5fb35b5a369a": 'biobarcoding/inputs_schema/clustalw_formly.json',
    "c87f58b6-cb06-4d39-a0b3-72c2705c5ae1": 'biobarcoding/inputs_schema/paup_parsimony_formly.json',
    "3e0240e8-b978-48a2-8fdd-9f31f4264064": 'biobarcoding/inputs_schema/pda_formly.json',
    "903a73a9-5a4e-4cec-b8fa-4fc9bd5ffab5": 'biobarcoding/inputs_schema/mafft_formly.json',
    "985c01ca-d9d2-4df5-a8b9-8a6da251d7d4": 'biobarcoding/inputs_schema/blast_formly.json',
    "ea647c4e-2063-4246-bd9a-42f6a57fb9ea": 'biobarcoding/inputs_schema/mrbayes_formly.json',
    "c55280d0-f916-4401-a1a4-bb26d8179fd7": 'biobarcoding/inputs_schema/clustalw+paup_parsimony_formly.json',
    "ce018826-7b20-4b70-b9b3-168c0ba46eec": 'biobarcoding/inputs_schema/paup_parsimony+pda_formly.json',
    "5b315dc5-ad12-4214-bb6a-bf013f0e4b8c": 'biobarcoding/inputs_schema/clustalw+paup_parsimony+pda_formly.json',
    "848b46b0-8602-42a2-a3fd-1b9be728d729": 'biobarcoding/inputs_schema/beast_formly.json'
}


def drop_pg_database(sa_str, database_name):
    db_connection_string = sa_str
    data_engine = sqlalchemy.create_engine(db_connection_string, echo=False)
    conn = data_engine.connect()
    conn.execute("commit")
    try:
        conn.execute("drop database " + database_name)
    except:
        pass
    conn.close()
    data_engine.dispose()


def create_pg_database_engine(sa_str, database_name, recreate_db=False):
    if recreate_db:
        drop_pg_database(sa_str, database_name)
    db_connection_string = sa_str
    data_engine = sqlalchemy.create_engine(db_connection_string, echo=False)
    conn = data_engine.connect()
    conn.execute("commit")
    try:
        conn.execute("create database " + database_name)
    except:
        pass
    conn.close()
    data_engine.dispose()
    db_connection_string = sa_str + database_name
    return sqlalchemy.create_engine(db_connection_string, echo=False, poolclass=QueuePool, pool_size=0, max_overflow=-1)


def load_table(sf, clazz, d):
    """
    Insert pairs (key, value) into a relational table
    It loads a dictionary "d" containing keys and values, into a relational table associated to the class "clazz",
    using the session factory "sf"

    :param sf:
    :param clazz:
    :param d:
    :return:
    """
    session = sf()
    for k, v in d.items():
        i = session.query(clazz).filter(clazz.uuid == k).first()
        if not i:
            ins = clazz()
            ins.uuid = k
            ins.name = v
            session.add(ins)
        else:
            i.name = v
    session.commit()
    sf.remove()


def create_or_update_entity(session, clazz, attributes: List[str], t: Union[List, Tuple], id_attr="uuid", update=False):
    uuid_idx = attributes.index(id_attr)
    if t[uuid_idx] is not None:
        i = session.query(clazz).filter(clazz.uuid == t[uuid_idx]).first()
    else:
        i = None
    modify_attributes = update
    if not i:
        entity = clazz()
        add_to_session = True
        modify_attributes = True
    else:
        entity = i
        add_to_session = False

    if modify_attributes:
        for i, f in enumerate(t):
            v = f
            attr = attributes[i]
            if isinstance(attr, tuple):
                foreign_clazz = attr[0]
                foreign_filter_field = attr[1]
                foreign_refer_field = attr[2]
                attr = attr[3]
                i2 = session.query(foreign_clazz).filter(getattr(foreign_clazz, foreign_filter_field) == v).first()
                v = getattr(i2, foreign_refer_field)
            else:
                if isinstance(f, dict):
                    v = getattr(entity, attr)
                    if v is not None:
                        v = {**v, **f}
            setattr(entity, attr, v)

    if add_to_session:
        session.add(entity)

    return entity


def load_table_extended(sf, clazz, attributes: List[str], values: List[Tuple], update=False):
    """
    Insert records into a relational table
    Loads records in "values" using "attributes" names, into a relational table associated to the class "clazz",
    using the session factory "sf".

    NOTE: "uuid" field MUST be present in "attributes"

    :param sf:
    :param clazz:
    :param attributes: a list of attribute names
    :param update: True if update of attributes is wanted
    :return:
    """
    session = sf()
    for t in values:
        create_or_update_entity(session, clazz, attributes, t, update=update)
    session.commit()
    sf.remove()


def load_many_to_many_table(sf, clazz, lclazz, rclazz, attributes: List[str], values: List[Tuple]):
    """
    load_many_to_many_table(sf, RoleIdentity, Role, Identity, ["role_id", "identity_id"], [["admins", "test_user"]])

    :param sf:
    :param clazz:
    :param lclazz:
    :param rclazz:
    :param attributes:
    :param values:
    :return:
    """
    session = sf()
    for t in values:
        # Find id of left and right sides
        left = session.query(lclazz).filter(lclazz.name == t[0]).first()
        right = session.query(rclazz).filter(rclazz.name == t[1]).first()
        i = session.query(clazz).filter(
            and_(getattr(clazz, attributes[0]) == left.id, getattr(clazz, attributes[1]) == right.id)).first()
        if not i:
            ins = clazz()
            setattr(ins, attributes[0], left.id)
            setattr(ins, attributes[1], right.id)
            session.add(ins)
    session.commit()
    sf.remove()


def load_status_checkers(sf, app):
    session = sf()
    from urllib.parse import urlparse
    from ..services import get_or_create
    # connections from config
    from sqlalchemy.exc import IntegrityError
    for k, v in app.config.items():
        if isinstance(v, str) and '://' in v:
            u = urlparse(v)
            try:
                get_or_create(session, StatusChecker, name=k, type=u.scheme, url=u.geturl())
                session.commit()
            except IntegrityError:
                session.rollback()
                print('The StatusChecker %s may already exist.' % k)
    # connections for compute resources
    for crs in session.query(ComputeResource).all():
        try:
            get_or_create(session, StatusChecker, name=crs.name, type='compute-resource', resource_id=crs.id)
            session.commit()
        except IntegrityError:
            session.rollback()
            print('The StatusChecker %s may already exist.' % crs.name)
    sf.remove()


def load_computing_resources(sf):
    session = sf()
    with open(get_global_configuration_variable("RESOURCES_CONFIG_FILE_PATH")) as json_file:
        resources_dict = json.load(json_file)
    for resource_uuid, resource_dict in resources_dict.items():
        jm_type = session.query(JobManagementType).filter(
            JobManagementType.name == resource_dict["job_management_type"]).first()
        r = session.query(ComputeResource).filter(ComputeResource.uuid == resource_uuid).first()
        if not r:
            r = ComputeResource()
            r.uuid = resource_uuid
            r.name = resource_dict["name"]
            session.add(r)
        # Overwrite
        r.jm_type = jm_type
        r.jm_location = resource_dict["job_management_location"]
        r.jm_credentials = resource_dict["job_management_credentials"]
        r.jm_params = resource_dict.get("job_management_params", {})

    session.commit()
    sf.remove()


def load_processes_in_computing_resources(sf):
    session = sf()
    with open(get_global_configuration_variable("RESOURCES_CONFIG_FILE_PATH")) as json_file:
        resources_dict = json.load(json_file)
    from ..rest import tm_processes
    for resource_uuid, resource_dict in resources_dict.items():
        for k, v in tm_processes.items():
            if v in RESOURCE_PROCESSES_DICT[resource_dict["name"]]:
                process = session.query(Process).filter(Process.uuid == k).first()
                resource = session.query(ComputeResource).filter(
                    ComputeResource.uuid == resource_uuid).first()
                r = session.query(ProcessInComputeResource).filter(
                    and_(ProcessInComputeResource.process_id == process.id,
                         ProcessInComputeResource.resource_id == resource.id)).first()
                if not r:
                    r = ProcessInComputeResource()
                    # r.uuid = local_uuid
                    r.native_process_id = v
                    r.process = process
                    r.resource = resource
                    session.add(r)
    session.commit()
    sf.remove()


def load_process_input_schema(sf):
    session = sf()
    for k, v in PROCESSES_INPUTS.items():
        process = session.query(Process).filter(Process.uuid == k).first()
        if process:
            if not process.schema_inputs:
                path = os.path.join(ROOT, v)
                with open(path, 'r') as f:
                    inputs = json.load(f)
                    process.schema_inputs = inputs
    session.commit()
    sf.remove()


def is_testing_enabled(flask_app):
    if "TESTING" in flask_app.config:
        if isinstance(flask_app.config["TESTING"], bool):
            testing = flask_app.config["TESTING"]
        else:
            testing = flask_app.config["TESTING"].lower() in ["true", "1"]
    else:
        testing = False
    return testing


def reset_database(flask_app):
    """
    Empty ALL data in the database !!!!

    Used in testing web services

    :return:
    """
    if is_testing_enabled(flask_app):
        connection2 = engine.connect()
        tables = ORMBase.metadata.tables
        table_existence = [engine.dialect.has_table(connection2, tables[t].name) for t in tables]
        connection2.close()
        if False in table_existence:
            ORMBase.metadata.bind = engine
            ORMBase.metadata.create_all()

    for tbl in reversed(ORMBase.metadata.sorted_tables):
        engine.execute(tbl.delete())


# #####################################################################################################################
# >>>> IDENTIFIERS CASE MANAGEMENT FUNCTIONS <<<<
# #####################################################################################################################

class CaseInsensitiveDict(collections.MutableMapping):
    """
    A dictionary with case insensitive Keys.
    Prepared also to support TUPLES as keys, required because compound keys are required
    """

    def __init__(self, data=None, **kwargs):
        from collections import OrderedDict
        self._store = OrderedDict()
        if data is None:
            data = {}
        self.update(data, **kwargs)

    def encode(self):
        return self.get_data()

    def get_original_data(self):
        return {casedkey: mappedvalue for casedkey, mappedvalue in self._store.values()}

    def get_data(self):
        return {key: self._store[key][1] for key in self._store}

    def __setitem__(self, key, value):
        # Use the lowercased key for lookups, but store the actual
        # key alongside the value.
        if not isinstance(key, tuple):
            self._store[key.lower()] = (key, value)
        else:
            self._store[tuple([k.lower() for k in key])] = (key, value)

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            return self._store[key.lower()][1]
        else:
            return self._store[tuple([k.lower() for k in key])][1]

    def __delitem__(self, key):
        if not isinstance(key, tuple):
            del self._store[key.lower()]
        else:
            del self._store[tuple([k.lower() for k in key])]

    def __iter__(self):
        return (casedkey for casedkey, mappedvalue in self._store.values())

    def __len__(self):
        return len(self._store)

    def lower_items(self):
        """Like iteritems(), but with all lowercase keys."""
        return (
            (lowerkey, keyval[1])
            for (lowerkey, keyval)
            in self._store.items()
        )

    def __contains__(self, key):  # "in" operator to check if the key is present in the dictionary
        if not isinstance(key, tuple):
            return key.lower() in self._store
        else:
            return tuple([k.lower() for k in key]) in self._store

    def __eq__(self, other):
        if isinstance(other, collections.Mapping):
            other = CaseInsensitiveDict(other)
        else:
            return NotImplemented
        # Compare insensitively
        return dict(self.lower_items()) == dict(other.lower_items())

    # Copy is required
    def copy(self):
        return CaseInsensitiveDict(self._store.values())

    def __repr__(self):
        return str(dict(self.items()))


def create_dictionary(case_sens=case_sensitive, multi_dict=False, data=dict()):
    """
    Factory to create dictionaries

    :param case_sens: True to create a case sensitive dictionary, False to create a case insensitive one
    :param multi_dict: True to create a "MultiDict", capable of storing several values
    :param data: Dictionary with which the new dictionary is initialized
    :return:
    """

    if not multi_dict:
        if case_sens:
            tmp = {}
            tmp.update(data)
            return tmp  # Normal, "native" dictionary
        else:
            return CaseInsensitiveDict(data)
    else:
        if case_sens:
            return MultiDict(data)
        else:
            return CIMultiDict(data)


def strcmp(s1, s2):
    """
    Compare two strings for equality or not, considering a flag for case sensitiveness or not

    It also removes leading and trailing whitespace from both strings, so it is not sensitive to this possible
    difference, which can be a source of problems

    :param s1:
    :param s2:
    :return:
    """
    # Handling empty or None strings
    if not s1:
        return True if not s2 else False
    if not s2:
        return False

    if case_sensitive:
        return s1.strip() == s2.strip()
    else:
        return s1.strip().lower() == s2.strip().lower()


# #####################################################################################################################
# >>>> PARTIAL RETRIEVAL DICTIONARY CLASS <<<<
# #####################################################################################################################

class PartialRetrievalDictionary:
    def __init__(self):
        # A dictionary of key-name to dictionaries, where the dictionaries are each of the values of the key and the
        # value is a set of IDs having that value
        # dict(key-name, dict(key-value, set(obj-IDs with that key-value))
        self._keys = {}
        # Dictionary from ID to the tuple (composite-key-elements dict, object)
        self._objs = {}
        self._rev_objs = {}  # From object to ID
        # Counter
        self._id_counter = 0

    def get(self, key, key_and_value=False, full_key=False, just_oid=False):
        """
        Retrieve one or more objects matching "key"
        If "key_and_value" is True, return not only the value, also matching key (useful for multiple matching keys)
        If "full_key" is True, zero or one objects should be the result
        :param key:
        :param full_key:
        :return: A list of matching elements
        """
        if True:
            # Lower case values
            # Keys can be all lower case, because they will be internal Key components, not specified by users
            if case_sensitive:
                key2 = {k.lower(): v for k, v in key.items()}
            else:
                key2 = {k.lower(): v if k.startswith("__") else v.lower() for k, v in key.items()}
        else:
            key2 = key

        sets = [self._keys.get(k, {}).get(v, set()) for k, v in key2.items()]

        # Find shorter set and Remove it from the list
        min_len = 1e30
        min_len_set_idx = None
        for i, s in enumerate(sets):
            if len(s) < min_len:
                min_len = len(s)
                min_len_set_idx = i
        min_len_set = sets[min_len_set_idx]
        del sets[min_len_set_idx]
        # Compute intersections
        result = min_len_set.intersection(*sets)
        if just_oid:
            return result

        # Obtain list of results
        if full_key and len(result) > 1:
            raise Exception("Zero or one results were expected. " + str(len(result) + " obtained."))
        if not key_and_value:
            return [self._objs[oid][1] for oid in result]
        else:
            return [self._objs[oid] for oid in result]

    def get_one(self, key, key_and_value=False, full_key=False, just_oid=False):
        results = self.get(key, key_and_value, full_key, just_oid)
        return results[0] if results else None

    def put(self, key, value):
        """
        Insert implies the key does not exist
        Update implies the key exists
        Upsert does not care

        :param key:
        :param value:
        :return:
        """
        ptype = 'i'  # 'i', 'u', 'ups' (Insert, Update, Upsert)
        if True:
            # Lower case values
            # Keys can be all lower case, because they will be internal Key components, not specified by users
            if case_sensitive:
                key2 = {k.lower(): v for k, v in key.items()}
            else:
                key2 = {k.lower(): v if k.startswith("__") else v.lower() if isinstance(v, str) else v for k, v in
                        key.items()}
        else:
            key2 = key
        # Arrays containing key: values "not-present" and "present"
        not_present = []  # List of tuples (dictionary of key-values, value to be stored)
        present = []  # List of sets storing IDs having same key-value
        for k, v in key2.items():
            d = self._keys.get(k, {})
            if len(d) == 0:
                self._keys[k] = d
            if v not in d:
                not_present.append((d, v))
            else:
                present.append(d.get(v))

        if len(not_present) > 0:
            is_new = True
        else:
            if len(present) > 1:
                is_new = len(present[0].intersection(*present[1:])) == 0
            elif len(present) == 1:
                is_new = len(present) == 0
            else:
                is_new = False

        # Insert, Update or Upsert
        if is_new:  # It seems to be an insert
            # Check
            if ptype == 'u':
                raise Exception("Key does not exist")
            # Insert
            if value in self._rev_objs:
                oid = self._rev_objs[value]
            else:
                self._id_counter += 1
                oid = self._id_counter
                self._objs[oid] = (key, value)
                self._rev_objs[value] = oid

            # Insert
            for d, v in not_present:
                s = set()
                d[v] = s
                s.add(oid)
            for s in present:
                s.add(oid)
        else:
            if ptype == 'i':
                raise Exception("Key '+" + str(key2) + "' already exists")
            # Update
            # Find the ID for the key
            res = self.get(key, just_oid=True)
            if len(res) != 1:
                raise Exception("Only one result expected")
            # Update value (key is the same, ID is the same)
            self._objs[res[0]] = value

    def delete(self, key):
        def delete_single(key):
            if True:
                # Lower case values
                # Keys can be all lower case, because they will be internal Key components, not specified by users
                if case_sensitive:
                    key2 = {k.lower(): v for k, v in key.items()}
                else:
                    key2 = {k.lower(): v if k.startswith("__") else v.lower() for k, v in key.items()}
            else:
                key2 = key

            # Get IDs
            oids = self.get(key, just_oid=True)
            if len(oids) > 0:
                # From key_i: value_i remove IDs (set difference)
                for k, v in key2.items():
                    d = self._keys.get(k, None)
                    if d:
                        s = d.get(v, None)
                        if s:
                            s2 = s.difference(oids)
                            d[v] = s2
                            if not s2:
                                del d[v]  # Remove the value for the key

                # Delete oids
                for oid in oids:
                    del self._objs[oid]

                return len(oids)
            else:
                return 0

        if isinstance(key, list):
            res_ = 0
            for k in key:
                res_ += delete_single(k)
            return res_
        else:
            return delete_single(key)

    def to_pickable(self):
        # Convert to a jsonpickable structure
        return dict(keys=self._keys, objs=self._objs, cont=self._id_counter)

    def from_pickable(self, inp):
        self._keys = inp["keys"]
        self._objs = {int(k): v for k, v in inp["objs"].items()}
        self._rev_objs = {v[1]: k for k, v in self._objs.items()}
        self._id_counter = inp["cont"]

        return self  # Allows the following: prd = PartialRetrievalDictionary().from_pickable(inp)
