from ..main import BasicService, get_orm


##
# ANNOTATION FORM SERVICE
##

class FormRelationshipService(BasicService):

	def __init__(self):
		super(FormRelationshipService, self).__init__()
		self.orm = get_orm('form_relationship')


##
# ANNOTATION FORM SERVICE
##

class RelationshipService(BasicService):

	def __init__(self):
		super(RelationshipService, self).__init__()
		self.orm = get_orm('relationship')
