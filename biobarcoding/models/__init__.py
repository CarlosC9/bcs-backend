from copy import deepcopy

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker, class_mapper, ColumnProperty, RelationshipProperty

DBSession = scoped_session(sessionmaker())


class BaseMixin(object):
    # query = DBSession.query_property()
    # @declared_attr
    # def __tablename__(cls):
    # return cls.__name__.lower()

    # def __new__(cls, *args, **kwargs):
    #     obj = super().__new__(cls)
    #     if '_sa_instance_state' not in obj.__dict__:
    #         obj._sa_instance_state = InstanceState(obj, obj._sa_class_manager)
    #     return obj

    def __deepcopy__(self, memo):
        cls = self.__class__
        print(str(cls))
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k in class_mapper(cls).iterate_properties:
            if isinstance(k, ColumnProperty):
                name = k.columns[0].name
                getattr(self, name, None)
            elif isinstance(k, RelationshipProperty):
                #if k.back_populates:
                name = k.strategy.key
                getattr(self, name, None)

        for k, v in self.__dict__.items():
            deepcopy(v, memo)
            # setattr(result, k, deepcopy(v, memo))
        return result


ORMBase = declarative_base(cls=BaseMixin)
