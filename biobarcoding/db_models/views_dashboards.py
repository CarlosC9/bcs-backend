import uuid as uuid
from sqlalchemy import Column, Integer, String, ForeignKey, BigInteger
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship, backref

from . import ORMBase, GUID
from .core import FunctionalObject, data_object_type_id

prefix = "vis_"


class View(FunctionalObject):
    __tablename__ = f"{prefix}views"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['view'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)
    name = Column(String, nullable=False)
    description = Column(String)
    type = Column(String(50))  # Map, Chart, PivotTable, etc.
    data = Column(JSONB)  # Maps: layers, visible attribute, opacity, color palette, legend, projection, etc.


class Dashboard(FunctionalObject):
    __tablename__ = f"{prefix}dashboards"
    __mapper_args__ = {
        'polymorphic_identity': data_object_type_id['dashboard'],
    }
    id = Column(BigInteger, ForeignKey(FunctionalObject.id), primary_key=True)
    name = Column(String, nullable=False)
    description = Column(String)
    type = Column(String(50))  # Time centered (like Grafana), Only Maps, heterogeneous, etc.
    # Information for the context of the dashboard. Each view should store its position and extents
    data = Column(JSONB)


class DashboardView(ORMBase):
    __tablename__ = f"{prefix}dashboard_views"
    id = Column(Integer, primary_key=True, autoincrement=True)
    dashboard_id = Column(BigInteger, ForeignKey(Dashboard.id), primary_key=True)
    view_id = Column(BigInteger, ForeignKey(View.id), primary_key=True)
    dashboard = relationship(Dashboard, foreign_keys=[dashboard_id], backref=backref("views", cascade="all, delete-orphan"))
    view = relationship(View, foreign_keys=[view_id], backref=backref("dashboards", cascade="all, delete-orphan"))
    view_in_dashboard = Column(JSONB)


"""
Operations
* CRUD view (map), list views of a type,
* CRUD dashboard, list dashboards,  
  - lista de elementos. No sólo vistas
  - Posición de elementos, instanciar vista (puede ser anónima, creada en el dashboard, o creada en otro lugar especializado)
  - Navegación entre cuadros de mandos?
  - Drill down, drill up, etc.
  - Eventos entre vistas (vistas deberán tener cierta interface para poder ser interconectadas)
"""
