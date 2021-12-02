import uuid

from sqlalchemy import Integer, Column, String, ForeignKey, JSON, BigInteger, LargeBinary
from sqlalchemy.orm import relationship, backref

from . import ORMBase, GUID
from .core import FunctionalObject

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


class FunctionalObjectInFile(ORMBase):
    __tablename__ = f"{prefix}fos_in_files"
    file_id = Column(BigInteger, ForeignKey(File.id), nullable=False, primary_key=True)
    file = relationship(File, backref=backref("bos", cascade="all, delete-orphan"))
    fos_id = Column(BigInteger, ForeignKey(FunctionalObject.id), nullable=False, primary_key=True)
    fos = relationship(FunctionalObject, backref=backref("files", cascade="all, delete-orphan"))
