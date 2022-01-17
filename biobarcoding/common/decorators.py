import functools
import gzip
import warnings
from typing import IO, Dict

from celery import Celery
from flask import after_this_request, request


# Some decorators from "https://wiki.python.org/moin/PythonDecoratorLibrary"


def rest_function(func):
    """
    Aggregation of decorators used in every REST function

    :param func:
    :return:
    """
    # Open-Close SQLAlchemy session
    # Control authentication
    # Flask-Session
    # Log
    #


def singleton(cls):
    """ Use class as singleton.

 Sample use:


@singleton
class Foo:
    def __new__(cls):
        cls.x = 10
        return object.__new__(cls)

    def __init__(self):
        assert self.x == 10
        self.x = 15

assert Foo().x == 15
Foo().x = 20
assert Foo().x == 20

    """

    cls.__new_original__ = cls.__new__

    @functools.wraps(cls.__new__)
    def singleton_new(cls, *args, **kw):
        it = cls.__dict__.get('__it__')
        if it is not None:
            return it

        cls.__it__ = it = cls.__new_original__(cls, *args, **kw)
        it.__init_original__(*args, **kw)
        return it

    cls.__new__ = singleton_new
    cls.__init_original__ = cls.__init__
    cls.__init__ = object.__init__

    return cls


def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.

=== Examples of use ===

@deprecated
def some_old_function(x,y):
    return x + y

class SomeClass:
    @deprecated
    def some_old_method(self, x,y):
        return x + y

    """

    def new_func(*args, **kwargs):
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)

    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func


class countcalls(object):
    """ Decorator that keeps track of the number of times a function is called. """

    __instances = {}

    def __init__(self, f):
        self.__f = f
        self.__numcalls = 0
        countcalls.__instances[f] = self

    def __call__(self, *args, **kwargs):
        self.__numcalls += 1
        return self.__f(*args, **kwargs)

    @staticmethod
    def count(f):
        """ Return the number of times the function f was called """
        return countcalls.__instances[f].__numcalls

    @staticmethod
    def counts():
        """ Return a dict of {function: # of calls} for all registered functions """
        return dict([(f, countcalls.count(f)) for f in countcalls.__instances])


class Memoize:
    """
    Cache of function calls (non-persistent, non-refreshable)
    """

    def __init__(self, fn):
        self.fn = fn
        self.memo = {}

    def __call__(self, *args):
        if args not in self.memo:
            self.memo[args] = self.fn(*args)
        return self.memo[args]


class Memoize2(object):
    """cache the return value of a method

    This class is meant to be used as a decorator of methods. The return value
    from a given method invocation will be cached on the instance whose method
    was invoked. All arguments passed to a method decorated with memoize must
    be hashable.

    If a memoized method is invoked directly on its class the result will not
    be cached. Instead the method will be invoked like a static method:
    class Obj(object):
        @memoize
        def add_to(self, arg):
            return self + arg
    Obj.add_to(1) # not enough arguments
    Obj.add_to(1, 2) # returns 3, result is not cached
    """

    def __init__(self, func):
        self.func = func

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self.func
        return functools.partial(self, obj)

    def __call__(self, *args, **kw):
        obj = args[0]
        try:
            cache = obj.__cache
        except AttributeError:
            cache = obj.__cache = {}
        key = (self.func, args[1:], frozenset(kw.items()))
        try:
            res = cache[key]
        except KeyError:
            res = cache[key] = self.func(*args, **kw)
        return res


def gzipped(f):
    """
    Decorator to ZIP the response in a RESTful call
    """

    @functools.wraps(f)
    def view_func(*args, **kwargs):
        @after_this_request
        def zipper(response):
            accept_encoding = request.headers.get('Accept-Encoding', '')

            if 'gzip' not in accept_encoding.lower():
                return response

            response.direct_passthrough = False

            if (response.status_code < 200 or
                    response.status_code >= 300 or
                    'Content-Encoding' in response.headers):
                return response
            gzip_buffer = IO()
            gzip_file = gzip.GzipFile(mode='wb',
                                      fileobj=gzip_buffer)
            gzip_file.write(response.data)
            gzip_file.close()

            response.data = gzip_buffer.getvalue()
            response.headers['Content-Encoding'] = 'gzip'
            response.headers['Vary'] = 'Accept-Encoding'
            response.headers['Content-Length'] = len(response.data)

            return response

        return f(*args, **kwargs)

    return view_func


class celery_wf(object):
    """
    A decorator to ease the definition of workflows where the next task depends on an outcome and the same task may
    repeat many times.
    It expects: a Celery instance, a Workflow definition, a Celery task name:

    @celery_wf(celery_app, wf1, "export")

    The definition (wf1 in the example) is a dictionary of task names (terminal tasks do not need to appear),
     where each entry contains a dictionary with possible jumps.
     A jump to itself is always allowed and is implicitly defined.

    Depending on returned values, the next task can be:
    * return None (or no "return" statement) -> Next task using the same input parameter.
      The next task is defined as "": "<next task name" in the workflow definition
    * return <ret> -> Go to next task passing <ret> as input
    * return <int:countdown>, <ret> -> Repeat the task, waiting <countdown> seconds and passing <ret> as parameter
    * return None, <ret> -> Repeat the task as soon as possible (no countdown), passing <ret> as parameter
    * return "<outcome in the workflow definition", <ret> -> Go to task matching the <outcome label> passing <ret> as input

    See "definitions.py" for examples of use.

    """

    def __init__(self, cel_app: Celery, wf_def: Dict, name: str):
        self.celery_app = cel_app
        self.wf_def = wf_def
        self.name = name

    @staticmethod
    def append_text(s: str):
        pass
        # file = "/home/rnebot/Downloads/borrame/log.txt"
        # with open(file, "a+") as f:
        #     f.write(f"{s}\n")

    def __call__(self, f):
        def wrapped_f(*args):
            # Before call
            celery_wf.append_text("DECOR - before task")
            res = f(*args)
            if isinstance(res, tuple):
                celery_wf.append_text("DECOR - returned tuple")
                if len(res) == 2:
                    result, task_ctx = res[0], res[1]
                else:
                    raise Exception(f"Tuple with {len(res)} elements returned: {res}")
            elif res is None:
                celery_wf.append_text("DECOR - returned None")
                result = ""
                task_ctx = args[0]
            else:
                celery_wf.append_text("DECOR - returned a single result")
                result = ""
                task_ctx = res

            # After call, program a new task
            countdown = 0
            if result is None or isinstance(result, int):
                if isinstance(result, int):
                    countdown = result
                    celery_wf.append_text(f"DECOR - reenter applying a countdown of {countdown}")
                next_task = self.name
            else:
                d = self.wf_def.get(self.name, {})
                next_task = d.get(result, None)
            if next_task:
                if countdown > 0:
                    self.celery_app.signature(next_task).apply_async(args=(task_ctx,), countdown=countdown)
                else:
                    self.celery_app.signature(next_task).apply_async(args=(task_ctx,))
            else:
                celery_wf.append_text(f"DECOR - no next task. FINISHED")

        return wrapped_f
