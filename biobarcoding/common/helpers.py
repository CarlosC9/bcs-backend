import sys

import jsonpickle


def serialize_from_object(obj):
    tmp = sys.getrecursionlimit()
    sys.setrecursionlimit(10000)
    tmp_str = jsonpickle.encode(obj)  # .encode("ascii")
    sys.setrecursionlimit(tmp)
    return tmp_str


def deserialize_to_object(s):
    tmp = sys.getrecursionlimit()
    sys.setrecursionlimit(10000)
    tmp_str = jsonpickle.decode(s)
    sys.setrecursionlimit(tmp)
    return tmp_str
