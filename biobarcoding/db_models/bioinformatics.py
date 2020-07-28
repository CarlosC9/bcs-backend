# BIOINFORMATIC data
import sqlalchemy as sa
from biobarcoding.db_models import ORMBase, GUID


class BioinformaticObject(ORMBase):
    __versioned__ = {}
    __tablename__ = 'bos'
    id = sa.Column(sa.BigInteger, primary_key=True, autoincrement=True)
    uuid = sa.Column(GUID, unique=True)
    name = sa.Column(sa.Unicode)
    content = sa.Column(sa.UnicodeText)


class Sequence(BioinformaticObject):
    pass


class MultipleSequenceAlignment(BioinformaticObject):
    pass


class PhylogeneticTree(BioinformaticObject):
    pass
