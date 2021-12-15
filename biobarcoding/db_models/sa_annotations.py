from sqlalchemy import Column, ForeignKey, UniqueConstraint, Boolean, Integer, BigInteger, String, Sequence
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship, backref

from . import ORMBase, GUID, ObjectType
from .core import FunctionalObject

# ANNOTATION FORM

prefix = "sa_annotation_"


class AnnotationFormItem(ORMBase):
    __tablename__ = f"{prefix}form_item"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    name = Column(String(80), nullable=False, unique=True)
    description = Column(String(512))
    cvterm_id = Column(Integer, unique=True)
    dbxref_id = Column(Integer)


class AnnotationFormItemObjectType(ORMBase):
    __tablename__ = f"{prefix}form_item_object_type"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    form_item_id = Column(Integer, ForeignKey(AnnotationFormItem.id), nullable=False)
    object_type_id = Column(Integer, ForeignKey(ObjectType.id), nullable=False)
    form_item = relationship(AnnotationFormItem, backref=backref("annotation_form_items", cascade="all, delete-orphan"))
    object_type = relationship(ObjectType, backref=backref("object_types", cascade="all, delete-orphan"))

    __table_args__ = (
        UniqueConstraint(form_item_id, object_type_id, name=__tablename__ + '_c1'),
    )


class AnnotationFormField(AnnotationFormItem):
    __tablename__ = f"{prefix}form_field"

    id = Column(BigInteger, ForeignKey(AnnotationFormItem.id, ondelete="CASCADE"), primary_key=True)
    type = Column(String(32), nullable=False)   # tag, attribute, relationship
    range = Column(JSONB)
    view_type = Column(String(32))  # check, radio, date, etc.
    multiple = Column(Boolean, default=False)

    __mapper_args__ = {
        'polymorphic_identity': 'annotation_form_field',
        'polymorphic_on': type
    }


class AnnotationFormTag(AnnotationFormField):
    __tablename__ = f"{prefix}form_tag"
    __mapper_args__ = {
        'polymorphic_identity': 'tag',
    }
    id = Column(BigInteger, ForeignKey(AnnotationFormField.id, ondelete="CASCADE"), primary_key=True)


class AnnotationFormAttribute(AnnotationFormField):
    __tablename__ = f"{prefix}form_attribute"
    __mapper_args__ = {
        'polymorphic_identity': 'attribute',
    }
    id = Column(BigInteger, ForeignKey(AnnotationFormField.id, ondelete="CASCADE"), primary_key=True)


class AnnotationFormRelationship(AnnotationFormField):
    __tablename__ = f"{prefix}form_relationship"
    __mapper_args__ = {
        'polymorphic_identity': 'relationship',
    }
    id = Column(BigInteger, ForeignKey(AnnotationFormField.id, ondelete="CASCADE"), primary_key=True)


class AnnotationFormTemplate(AnnotationFormItem):
    __tablename__ = f"{prefix}form_template"

    id = Column(BigInteger, ForeignKey(AnnotationFormItem.id, ondelete="CASCADE"), primary_key=True)


class AnnotationFormTemplateField(ORMBase):
    __tablename__ = f"{prefix}form_template_field"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    form_field_id = Column(Integer, ForeignKey(AnnotationFormField.id), nullable=False)
    form_template_id = Column(Integer, ForeignKey(AnnotationFormTemplate.id), nullable=False)
    form_field = relationship(AnnotationFormField, backref=backref("annotation_form_fields", cascade="all, delete-orphan"))
    form_template = relationship(AnnotationFormTemplate, backref=backref("annotation_form_templates", cascade="all, delete-orphan"))
    name = Column(String(80))
    rank = Column(Integer, Sequence('annotation_form_rank_seq'), nullable=False)
    # rank = Column(Integer, nullable=False, default=select([func.max(1, func.max(rank))]))

    __table_args__ = (
        UniqueConstraint(form_template_id, rank, name=__tablename__ + '_c1'),
    )


# ANNOTATION INSTANCE

class AnnotationItem(ORMBase):
    __tablename__ = f"{prefix}item"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    object_uuid = Column(GUID, ForeignKey(FunctionalObject.uuid, ondelete="CASCADE"), nullable=False)
    type = Column(String(80), nullable=False)    # text, form, field
    name = Column(String(80), nullable=False)
    rank = Column(Integer, Sequence('annotation_rank_seq'), nullable=False)
    # rank = Column(Integer, nullable=False, default=select([func.max(1, func.max(rank))]))

    __table_args__ = (
        UniqueConstraint(object_uuid, rank, name=__tablename__ + '_c1'),
        UniqueConstraint(object_uuid, name, name=__tablename__ + '_c2'),
    )
    __mapper_args__ = {
        'polymorphic_identity': 'annotation_item',
        'polymorphic_on': type
    }


class AnnotationText(AnnotationItem):
    __tablename__ = f"{prefix}text"
    __mapper_args__ = {
        'polymorphic_identity': 'text',
    }
    id = Column(BigInteger, ForeignKey(AnnotationItem.id, ondelete="CASCADE"), primary_key=True)
    value = Column(String(512))


class AnnotationTemplate(AnnotationItem):
    __tablename__ = f"{prefix}template"
    __mapper_args__ = {
        'polymorphic_identity': 'template',
    }
    id = Column(BigInteger, ForeignKey(AnnotationItem.id, ondelete="CASCADE"), primary_key=True)
    form_template_id = Column(Integer, ForeignKey(AnnotationFormTemplate.id), nullable=False)
    form_template = relationship(AnnotationFormTemplate, backref=backref("annotation_templates", cascade="all, delete-orphan"))
    value = Column(JSONB)


class AnnotationField(AnnotationItem):
    __tablename__ = f"{prefix}field"
    __mapper_args__ = {
        'polymorphic_identity': 'field',
    }
    id = Column(BigInteger, ForeignKey(AnnotationItem.id, ondelete="CASCADE"), primary_key=True)
    form_field_id = Column(Integer, ForeignKey(AnnotationFormField.id), nullable=False)
    form_field = relationship(AnnotationFormField, backref=backref("annotation_fields", cascade="all, delete-orphan"))
    value = Column(JSONB)
