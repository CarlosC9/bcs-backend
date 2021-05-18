from biobarcoding.rest import Issue, IType


def create(**kwargs):
    issues, status = [Issue(IType.INFO, f'CREATE views: dummy successfully completed.')], 200
    return issues, None, 0, status


def read(view=None, **kwargs):
    issues, status = [Issue(IType.INFO, f'READ views({view}): dummy successfully completed.')], 200
    return issues, None, 0, status


def update(view, **kwargs):
    issues, status = [Issue(IType.INFO, f'UPDATE views({view}): dummy successfully completed.')], 200
    return issues, None, status


def delete(view=None, **kwargs):
    issues, status = [Issue(IType.INFO, f'DELETE views({view}): dummy successfully completed.')], 200
    return issues, None, status
