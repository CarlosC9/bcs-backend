import collections
import json
from pathlib import Path
from typing import Any, Dict
import uuid

import numpy as np
import pandas
from multidict import MultiDict, CIMultiDict


class Encodable:
    """
    Abstract class with the method encode() that should be implemented by a subclass to be encoded into JSON
    using the json.dumps() method together with the option cls=CustomEncoder.
    """

    def encode(self) -> Dict[str, Any]:
        raise NotImplementedError("users must define encode() to use this base class")

    @staticmethod
    def parents_encode(obj: "Encodable", cls: type) -> Dict[str, Any]:
        """
        Get the state of all "cls" parent classes for the selected instance "obj"
        :param obj: The instance. Use "self".
        :param cls: The base class which parents we want to get. Use "__class__".
        :return: A dictionary with the state of the instance "obj" for all inherited classes.

        """
        d = {}
        for parent in cls.__bases__:
            if issubclass(parent, Encodable) and parent is not Encodable:
                d.update(parent.encode(obj))
        return d


class CaseInsensitiveDict(collections.MutableMapping, Encodable):
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


def create_dictionary(case_sens=True, multi_dict=False, data=dict()):
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


def _json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""
    from datetime import datetime
    from uuid import UUID
    from ..db_models.jobs import StatusEnum

    if isinstance(obj, datetime):
        serial = obj.isoformat()
        return serial
    elif isinstance(obj, CaseInsensitiveDict):
        return str(obj)
    elif isinstance(obj, np.int64):
        return int(obj)
    elif isinstance(obj, pandas.DataFrame):
        return obj.to_dict("records")
    elif isinstance(obj, UUID):
        return str(obj)
    elif isinstance(obj, StatusEnum):
        return str(obj)
    elif hasattr(obj, "Schema"):
        return getattr(obj, "Schema")().dump(obj)  # "marshmallow"
    raise TypeError(f"Type {type(obj)} not serializable")


JSON_INDENT = 4
ENSURE_ASCII = False
ROOT = str(Path(__file__).parent.parent.parent)


def generate_json(o):
    return json.dumps(o,
                      default=_json_serial,
                      sort_keys=True,
                      indent=JSON_INDENT,
                      ensure_ascii=ENSURE_ASCII,
                      separators=(',', ': ')
                      ) if o else None


def check_pid_running(pid):
    import psutil
    try:
        proc = psutil.Process(pid)
        if proc.status() == psutil.STATUS_ZOMBIE:
            proc.kill()
            return False
    except psutil.NoSuchProcess:
        return False
    else:
        return True
