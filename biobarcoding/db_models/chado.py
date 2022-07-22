from sqlalchemy import Column, Integer, ForeignKey, BigInteger, Sequence, String
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.sql.functions import coalesce
from sqlalchemy_continuum.transaction import TransactionBase

from . import ORMBaseChado


##
# EXTENDED
##

class Analysis(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "analysis"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class AnalysisRelationship(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "analysis_relationship"
    __table_args__ = {'extend_existing': True, 'autoload': True}

    subject_id = Column(Integer, ForeignKey(Analysis.analysis_id))
    subject = relationship(Analysis, foreign_keys=[subject_id],
                           backref=backref("related_objects", cascade="all, delete-orphan"))
    object_id = Column(Integer, ForeignKey(Analysis.analysis_id))
    object = relationship(Analysis, foreign_keys=[object_id],
                          backref=backref("related_subjects", cascade="all, delete-orphan"))


class Organism(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "organism"
    __table_args__ = {'extend_existing': True, 'autoload': True}

    @hybrid_property
    def name(self):
        return self.genus + ' ' + self.species + coalesce(' ' + self.infraspecific_name, '')


class Phylotree(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "phylotree"
    __table_args__ = {'extend_existing': True, 'autoload': True}

    analysis = relationship(Analysis, backref=backref("phylotrees", cascade="all, delete-orphan"))


class AnalysisFeature(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "analysisfeature"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class Cv(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "cv"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class Cvterm(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "cvterm"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class Db(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "db"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class Dbxref(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "dbxref"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class Feature(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "feature"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class Featureloc(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "featureloc"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class Phylonode(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "phylonode"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class PhylonodeOrganism(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "phylonode_organism"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class Phylonodeprop(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "phylonodeprop"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class Stock(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "stock"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class StockFeature(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "stock_feature"
    __table_args__ = {'extend_existing': True, 'autoload': True}


class Stockcollection(ORMBaseChado):
    __versioned__ = {}
    __tablename__ = "stockcollection"
    __table_args__ = {'extend_existing': True, 'autoload': True}


##
# AUTOLOADED
##

class AnalysisCvterm(ORMBaseChado):
    __tablename__ = "analysis_cvterm"
    __table_args__ = {'autoload': True}


class AnalysisDbxref(ORMBaseChado):
    __tablename__ = "analysis_dbxref"
    __table_args__ = {'autoload': True}


class Analysisprop(ORMBaseChado):
    __tablename__ = "analysisprop"
    __table_args__ = {'autoload': True}


class FeatureCvterm(ORMBaseChado):
    __tablename__ = "feature_cvterm"
    __table_args__ = {'autoload': True}


class Featureprop(ORMBaseChado):
    __tablename__ = "featureprop"
    __table_args__ = {'autoload': True}


class Phylotreeprop(ORMBaseChado):
    __tablename__ = "phylotreeprop"
    __table_args__ = {'autoload': True}


class StockCvterm(ORMBaseChado):
    __tablename__ = "stock_cvterm"
    __table_args__ = {'autoload': True}


class StockcollectionStock(ORMBaseChado):
    __tablename__ = "stockcollection_stock"
    __table_args__ = {'autoload': True}


class Stockcollectionprop(ORMBaseChado):
    __tablename__ = "stockcollectionprop"
    __table_args__ = {'autoload': True}


class Stockprop(ORMBaseChado):
    __tablename__ = "stockprop"
    __table_args__ = {'autoload': True}


class Transaction(ORMBaseChado, TransactionBase):
    """
    Credits to SQLAlchemy-Continuum
    @source: https://github.com/kvesteri/sqlalchemy-continuum/blob/master/sqlalchemy_continuum/transaction.py#L128
    """
    __tablename__ = 'transaction'

    id = Column(BigInteger, Sequence('transaction_id_seq'), primary_key=True, autoincrement=True)
    remote_addr = Column(String(50))

    def __repr__(self):
        fields = ['id', 'issued_at', 'user']
        from collections import OrderedDict
        field_values = OrderedDict(
            (field, getattr(self, field))
            for field in fields
            if hasattr(self, field)
        )
        import six
        return '<Transaction %s>' % ', '.join(
            (
                '%s=%r' % (field, value)
                if not isinstance(value, six.integer_types)
                # We want the following line to ensure that longs get
                # shown without the ugly L suffix on python 2.x
                # versions
                else '%s=%d' % (field, value)
                for field, value in field_values.items()
            )
        )