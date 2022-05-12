from . import SysService
from ..main import get_orm

##
# SUBSYSTEM STATUS TOOLS
##


class Service(SysService):

	def __init__(self):
		super(Service, self).__init__()
		self.orm = get_orm('status_checkers')

	def check_values(self, **values) -> dict:
		from ...services import secure_url
		values['url'] = secure_url(values.get('url', ''))
		return super(Service, self).check_values(**values)

	def after_create(self, new_object, **values):
		new_object.check_status()
		return super(Service, self).after_create(new_object, **values)

	def after_update(self, new_object, **values):
		new_object.check_status()
		return super(Service, self).after_update(new_object, **values)