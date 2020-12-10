import datetime

from sqlalchemy import Integer, Column, String, Text, ForeignKey, JSON, DateTime, BigInteger, LargeBinary
from sqlalchemy.orm import relationship, backref

from biobarcoding.db_models import ORMBase, GUID
import uuid

from biobarcoding.db_models.bioinformatics import BioinformaticObject

prefix = "fs_"


class FileSystemObject(ORMBase):
    """  A file system object is a File or a Folder"""
    __versioned__ = {}
    __tablename__ = f"{prefix}objects"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    uuid = Column(GUID, unique=True, default=uuid.uuid4)
    fso_type = Column(String(6), nullable=False)
    name = Column(String(1024))
    full_name = Column(String(2048))
    __mapper_args__ = {'polymorphic_on': fso_type}


class File(FileSystemObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}files"
    id = Column(BigInteger, ForeignKey(FileSystemObject.id), primary_key=True)
    content_type = Column(String(256))
    embedded_content = Column(LargeBinary)
    content_size = Column(Integer)
    content_location = Column(JSON)  # Location and maybe credentials
    __mapper_args__ = {
        'polymorphic_identity': 'file'
    }


class Folder(FileSystemObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}folders"
    id = Column(BigInteger, ForeignKey(FileSystemObject.id), primary_key=True)
    __mapper_args__ = {
        'polymorphic_identity': 'folder'
    }


FileSystemObject.parent_id = Column(BigInteger, ForeignKey(Folder.id), nullable=True, primary_key=False)
FileSystemObject.parent = relationship(Folder, foreign_keys=[FileSystemObject.parent_id],
                                       backref=backref("children", cascade="all, delete-orphan"))


class BioinformaticObjectInFile(ORMBase):
    __tablename__ = f"{prefix}bos_in_files"
    file_id = Column(BigInteger, ForeignKey(File.id), nullable=False, primary_key=True)
    file = relationship(File)
    bos_id = Column(BigInteger, ForeignKey(BioinformaticObject.id), nullable=False, primary_key=True)
    bos = relationship(BioinformaticObject, backref=backref("files", cascade="all, delete-orphan"))
