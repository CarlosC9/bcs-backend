import json
import numpy as np


def _json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""
    from datetime import datetime

    if isinstance(obj, datetime):
        serial = obj.isoformat()
        return serial
    elif isinstance(obj, np.int64):
        return int(obj)
    raise TypeError(f"Type {type(obj)} not serializable")

JSON_INDENT = 4
ENSURE_ASCII = False


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