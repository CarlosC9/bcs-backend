from . import celery_app

"""
Celery tasks CANNOT be debugged in Celery!! (the code is run in a separate process; of course they can be debugged "off-line")
"""

# Send messages
# Refresh queues of jobs (at different computing resources)


@celery_app.task
def add(x, y):
    print(f"NORMAL CELERY TASK: {x}+{y}={x + y}")
    return x + y


@celery_app.task
def periodic_sum(x, y):
    print(f"PERIODIC CELERY TASK: {x}+{y}={x + y}")
    return x + y
