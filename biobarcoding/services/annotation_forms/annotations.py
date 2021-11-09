from ..main import SimpleAuxService, get_orm


##
# ANNOTATION TOOLS
##

class AuxService(SimpleAuxService):
    # TODO:
    #  if an existent field allow multiple instances
    #  if the value for a field is valid by its range
    #  if the value for a field is valid by its type (tag, attribute, relationship)

    def __init__(self):
        super(AuxService, self).__init__()
        self.orm = get_orm('annotations')

    def read(self, **kwargs):
        content, count = self.get_query(**kwargs)
        if kwargs.get('id') or kwargs.get('object_uuid'):
            content = content.first()
        else:
            content = content.all()
        return content, count