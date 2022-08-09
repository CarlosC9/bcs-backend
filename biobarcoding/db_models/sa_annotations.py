from sqlalchemy import Column, ForeignKey, UniqueConstraint, Boolean, Integer, BigInteger, String, Sequence, Text
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship, validates

from . import ORMBase, GUID, ObjectType
from .core import FunctionalObject, data_object_type_id
from .files import File

prefix = "sa_annotation_"

data_object_type_id.update({
    "annotation_item": 20,
    "annotation_template": 21,
    "annotation_field": 22,
    "annotation_text": 23,
    # "relationship": 24,
})


# ANNOTATION FORM DEFINITION


class AnnotationFormItem(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}form_item"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    name = Column(String(80), nullable=False, unique=True)
    description = Column(String(1024))
    standard = Column(String(32))   # 'bibtex', 'darwin core', 'dublin core', ...
    cvterm_id = Column(Integer, unique=True)
    dbxref_id = Column(Integer)
    type = Column(String(32), nullable=False)   # template, field
    unique = Column(Boolean, nullable=False, default=False)   # TODO: validation
    object_types = relationship("AnnotationFormItemObjectType", back_populates="form_item", cascade="all, delete-orphan")

    __mapper_args__ = {
        'polymorphic_identity': 'annotation_form_item',
        'polymorphic_on': type
    }


class AnnotationFormItemObjectType(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}form_item_object_type"

    form_item_id = Column(BigInteger, ForeignKey(AnnotationFormItem.id, ondelete="CASCADE"), primary_key=True)
    object_type_id = Column(Integer, ForeignKey(ObjectType.id, ondelete="CASCADE"), primary_key=True)
    form_item = relationship(AnnotationFormItem, back_populates="object_types")


class AnnotationFormField(AnnotationFormItem):
    __versioned__ = {}
    __tablename__ = f"{prefix}form_field"
    __mapper_args__ = {
        'polymorphic_identity': 'field',
    }
    id = Column(BigInteger, ForeignKey(AnnotationFormItem.id, ondelete="CASCADE"), primary_key=True)
    range = Column(JSONB)
    view_type = Column(String(32), default='annotation-text')  # check, radio, date, etc.
    schema = Column(JSONB)  # Schema of the data defined by the relationship
    annotation_form_templates = relationship("AnnotationFormTemplateField", back_populates="form_field")
    annotation_items = relationship("AnnotationField", back_populates="form_field")


class AnnotationFormTemplate(AnnotationFormItem):
    __versioned__ = {}
    __tablename__ = f"{prefix}form_template"
    __mapper_args__ = {
        'polymorphic_identity': 'template',
    }
    id = Column(BigInteger, ForeignKey(AnnotationFormItem.id, ondelete="CASCADE"), primary_key=True)
    annotation_form_fields = relationship("AnnotationFormTemplateField", back_populates="form_template")
    annotation_items = relationship("AnnotationTemplate", back_populates="form_template")


class AnnotationFormTemplateField(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}form_template_field"

    id = Column(BigInteger, primary_key=True, unique=True, autoincrement=True)
    form_template_id = Column(BigInteger, ForeignKey(AnnotationFormTemplate.id, ondelete="CASCADE"), primary_key=True)
    form_field_id = Column(BigInteger, ForeignKey(AnnotationFormField.id, ondelete="CASCADE"), primary_key=True)
    form_template = relationship(AnnotationFormTemplate, back_populates="annotation_form_fields")
    form_field = relationship(AnnotationFormField, back_populates="annotation_form_templates")
    name = Column(String(80))
    rev_name = Column(String(80))
    rank = Column(Integer, Sequence('annotation_form_rank_seq'), primary_key=True)
    # rank = Column(Integer, nullable=False, default=select([func.max(1, func.max(rank))]))
    required = Column(Boolean, nullable=False, default=False)   # TODO: validation
    relationships = relationship("AnnotationRelationship",
                                 primaryjoin="AnnotationRelationship.form_rl_id == AnnotationFormTemplateField.id",
                                 back_populates="form_rl", cascade="all, delete-orphan")

    @validates('form_field')
    def validate_field(self, key, field):
        # doesn't work when assignment
        if field.unique:
            assert field.id not in [f.id for f in self.form_field]
        return field

    __table_args__ = (
        UniqueConstraint(form_template_id, rank, name=__tablename__ + '_c1'),
        # ForeignKeyConstraint([form_field_id, unique], [AnnotationFormItem.id, AnnotationFormItem.unique]),
        # Index(__tablename__ + '_c2',
        #       form_template_id, form_field_id,
        #       unique=True,
        #       postgresql_where=Column("unique=1")),
        # CheckConstraint(form_template.object_types.intersect(form_field.object_types), name=__tablename__ + '_c2'),
    )


# ANNOTATION INSTANCE


class AnnotationItem(FunctionalObject):
    __versioned__ = {}
    __tablename__ = f"{prefix}item"

    id = Column(BigInteger, ForeignKey(FunctionalObject.id, ondelete="CASCADE"), primary_key=True)
    type = Column(String(80), nullable=False)    # text, form, field
    # name = Column(String(80), nullable=False)
    file_id = Column(BigInteger, ForeignKey(File.id))
    file = relationship(File)
    objects = relationship("AnnotationItemFunctionalObject", back_populates="annotation", cascade="all, delete-orphan")

    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['annotation_item'],
    }


class AnnotationItemFunctionalObject(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}item_functional_object"

    object_uuid = Column(GUID, ForeignKey(FunctionalObject.uuid, ondelete="CASCADE"), primary_key=True)
    annotation_id = Column(Integer, ForeignKey(AnnotationItem.id, ondelete="CASCADE"), primary_key=True)
    object = relationship(FunctionalObject, foreign_keys=[object_uuid], backref="annotations")
    annotation = relationship(AnnotationItem, foreign_keys=[annotation_id], back_populates="objects")
    rank = Column(Integer, Sequence('annotation_rank_seq'), nullable=False)
    # rank = Column(Integer, nullable=False, default=select([func.max(1, func.max(rank))]))

    __table_args__ = (
        UniqueConstraint(annotation_id, object_uuid, name=__tablename__ + '_c1'),
        UniqueConstraint(object_uuid, rank, name=__tablename__ + '_c2'),
    )

    # @validates('annotation')
    # def validate_field(self, key, annotation):
    #     # doesn't work when assignment
    #     if hasattr(annotation, 'form_field') and annotation.form_field.unique:
    #         assert annotation.id not in [f.id for f in self.annotation]
    #     if hasattr(annotation, 'form_template') and annotation.form_template.unique:
    #         assert annotation.id not in [f.id for f in self.annotation]
    #     return annotation


class AnnotationText(AnnotationItem):
    __versioned__ = {}
    __tablename__ = f"{prefix}text"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['annotation_text'],
    }

    id = Column(BigInteger, ForeignKey(AnnotationItem.id, ondelete="CASCADE"), primary_key=True)
    value = Column(Text, unique=True)

    def __init__(self, **kwargs):
        super(AnnotationText, self).__init__(**kwargs)
        self.type = 'text'


class AnnotationTemplate(AnnotationItem):
    __versioned__ = {}
    __tablename__ = f"{prefix}template"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['annotation_template'],
    }

    id = Column(BigInteger, ForeignKey(AnnotationItem.id, ondelete="CASCADE"), primary_key=True)
    form_template_id = Column(BigInteger, ForeignKey(AnnotationFormTemplate.id, ondelete="CASCADE"), nullable=False)
    form_template = relationship(AnnotationFormTemplate, back_populates="annotation_items")
    value = Column(JSONB)

    __table_args__ = (
        UniqueConstraint(form_template_id, value, name=__tablename__ + '_c1'),
    )

    def __init__(self, **kwargs):
        super(AnnotationTemplate, self).__init__(**kwargs)
        self.type = 'template'


class AnnotationField(AnnotationItem):
    __versioned__ = {}
    __tablename__ = f"{prefix}field"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['annotation_field'],
    }

    id = Column(BigInteger, ForeignKey(AnnotationItem.id, ondelete="CASCADE"), primary_key=True)
    form_field_id = Column(BigInteger, ForeignKey(AnnotationFormField.id, ondelete="CASCADE"), nullable=False)
    form_field = relationship(AnnotationFormField, back_populates="annotation_items")
    value = Column(JSONB)

    __table_args__ = (
        UniqueConstraint(form_field_id, value, name=__tablename__ + '_c1'),
    )

    def __init__(self, **kwargs):
        super(AnnotationField, self).__init__(**kwargs)
        self.type = 'field'


class AnnotationRelationship(ORMBase):
    __versioned__ = {}
    __tablename__ = f"{prefix}relationship"

    id = Column(BigInteger, primary_key=True, autoincrement=True)
    subject_id = Column(BigInteger, ForeignKey(FunctionalObject.id, ondelete="CASCADE"), nullable=False)
    subject = relationship(FunctionalObject, foreign_keys=[subject_id], backref="related_objects")
    object_id = Column(BigInteger, ForeignKey(FunctionalObject.id, ondelete="CASCADE"), nullable=False)
    object = relationship(FunctionalObject, foreign_keys=[object_id], backref="related_subjects")
    form_rl_id = Column(BigInteger, ForeignKey(AnnotationFormTemplateField.id, ondelete="CASCADE"), nullable=False)
    form_rl = relationship(AnnotationFormTemplateField, foreign_keys=[form_rl_id])
        # , back_populates="relationships")
    form_rev_rl_id = Column(BigInteger, ForeignKey(AnnotationFormTemplateField.id))
    form_rev_rl = relationship(AnnotationFormTemplateField, foreign_keys=[form_rev_rl_id])
        # , back_populates="rev_relationships")
    value = Column(JSONB)

    __table_args__ = (
        UniqueConstraint(subject_id, object_id, form_rl_id, name=__tablename__ + '_c1'),
    )
