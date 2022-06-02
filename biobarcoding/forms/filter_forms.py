from ..db_models.chado import Cvterm
from ..db_models.core import CaseStudy, CProcess
from ..db_models.hierarchies import HierarchyNode, Hierarchy
from ..rest import h_subjects_name, h_sources_name, h_crs_name
from ..db_models import DBSession, DBSessionChado


def __getTypes(type):
    if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
        from ..db_models.chado import Feature as ORM
    elif type == 'phylotree' or type == 'phylotrees':
        from ..db_models.chado import Phylotree as ORM
    elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
        from ..db_models.chado import Stock as ORM
    else:
        return None
    ids = DBSessionChado.query(ORM.type_id).distinct(ORM.type_id).subquery()
    return DBSessionChado.query(Cvterm).filter(Cvterm.cvterm_id.in_(ids)).all()


def __getCvterms(type):
    if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
        from ..db_models.chado import FeatureCvterm as ORM
    elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
        from ..db_models.chado import AnalysisCvterm as ORM
    elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
        from ..db_models.chado import StockCvterm as ORM
    elif type == 'annotation_form_template' or type == 'annotation_form_templates' \
            or type == 'annotation_form_field' or type == 'annotation_form_fields':
        from ..db_models.sa_annotations import AnnotationFormItem as ORM
        ids = DBSession.query(ORM.cvterm_id).distinct(ORM.cvterm_id).all()
        return DBSessionChado.query(Cvterm).filter(Cvterm.cvterm_id.in_(ids)).all()
    else:
        return None
    ids = DBSessionChado.query(ORM.cvterm_id).distinct(ORM.cvterm_id).subquery()
    return DBSessionChado.query(Cvterm).filter(Cvterm.cvterm_id.in_(ids)).all()


def __getProps(type):
    if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
        from ..db_models.chado import Featureprop as ORM
    elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
        from ..db_models.chado import Analysisprop as ORM
    elif type == 'phylotree' or type == 'phylotrees':
        from ..db_models.chado import Phylotreeprop as ORM
    elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
        from ..db_models.chado import Stockprop as ORM
    else:
        return None
    ids = DBSessionChado.query(ORM.type_id).distinct(ORM.type_id).subquery()
    return DBSessionChado.query(Cvterm).filter(Cvterm.cvterm_id.in_(ids)).all()


def __getDbxref(type):
    from ..db_models.chado import Dbxref
    if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
        from ..db_models.chado import Feature as ORM
    elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
        from ..db_models.chado import AnalysisDbxref as ORM
    elif type == 'phylotree' or type == 'phylotrees':
        from ..db_models.chado import Phylotree as ORM
    elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
        from ..db_models.chado import Stock as ORM
    elif type == 'annotation_form_template' or type == 'annotation_form_templates' \
            or type == 'annotation_form_field' or type == 'annotation_form_fields':
        from ..db_models.sa_annotations import AnnotationFormItem as ORM
        ids = DBSession.query(ORM.dbxref_id).distinct(ORM.dbxref_id).all()
        return DBSessionChado.query(Dbxref).filter(Dbxref.dbxref_id.in_(ids)).all()
    else:
        return None
    ids = DBSessionChado.query(ORM.dbxref_id).distinct(ORM.dbxref_id).subquery()
    return DBSessionChado.query(Dbxref).filter(Dbxref.dbxref_id.in_(ids)).all()


def __getOrganisms(type):
    from ..db_models.chado import Feature as ORM
    if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
        ids = DBSessionChado.query(ORM.organism_id)
    elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
        # AnalysisFeature > Feature
        from ..db_models.chado import AnalysisFeature
        ids = DBSessionChado.query(AnalysisFeature.feature_id).subquery()
        ids = DBSessionChado.query(ORM.organism_id).filter(ORM.feature_id.in_(ids))
    elif type == 'phylotree' or type == 'phylotrees':
        # Phylonode > Feature
        from ..db_models.chado import Phylonode
        ids = DBSessionChado.query(Phylonode.feature_id).subquery()
        ids = DBSessionChado.query(ORM.organism_id).filter(ORM.feature_id.in_(ids))
    elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
        from ..db_models.chado import Stock as ORM
        ids = DBSessionChado.query(ORM.organism_id)
    elif type == 'taxonomy' or type == 'taxonomies':
        from ..db_models.chado import PhylonodeOrganism as ORM
        ids = DBSessionChado.query(ORM.organism_id)
    else:
        return None
    ids = ids.distinct(ORM.organism_id).all()
    from ..services.bio.meta.organisms import Service as OrgService
    return OrgService().read(filter={'organism_id': ids})[0]    # info attachment needed


def __getAnalyses(type):
    if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
        from ..db_models.chado import AnalysisFeature as ORM
    elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
        from ..db_models.chado import Analysis as ORM
    elif type == 'phylotree' or type == 'phylotrees':
        from ..db_models.chado import Phylotree as ORM
    else:
        return None
    ids = DBSessionChado.query(ORM.analysis_id).distinct(ORM.analysis_id).subquery()
    from ..db_models.chado import Analysis
    return DBSessionChado.query(Analysis).filter(Analysis.analysis_id.in_(ids)).all()


def __getAnnotationFormTemplates(type):
    if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
        from ..db_models.bioinformatics import Sequence as ORM
    elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
        from ..db_models.bioinformatics import MultipleSequenceAlignment as ORM
    elif type == 'phylotree' or type == 'phylotrees':
        from ..db_models.bioinformatics import PhylogeneticTree as ORM
    elif type == 'annotation_form_field' or type == 'annotation_form_fields':
        from ..db_models.sa_annotations import AnnotationFormTemplate, AnnotationFormTemplateField as ORM
        return DBSession.query(AnnotationFormTemplate).join(ORM).all()
    else:
        return None
    from ..db_models.sa_annotations import AnnotationFormTemplate, AnnotationTemplate, AnnotationItemFunctionalObject
    return DBSession.query(AnnotationFormTemplate).join(AnnotationTemplate).join(AnnotationItemFunctionalObject).join(ORM).all()


def __getAnnotationFormFields(type):
    if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
        from ..db_models.bioinformatics import Sequence as ORM
    elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
        from ..db_models.bioinformatics import MultipleSequenceAlignment as ORM
    elif type == 'phylotree' or type == 'phylotrees':
        from ..db_models.bioinformatics import PhylogeneticTree as ORM
    elif type == 'annotation_form_template' or type == 'annotation_form_templates':
        from ..db_models.sa_annotations import AnnotationFormField, AnnotationFormTemplateField as ORM
        return DBSession.query(AnnotationFormField).join(ORM).all()
    else:
        return None
    from ..db_models.sa_annotations import AnnotationFormField, AnnotationField, AnnotationItemFunctionalObject
    return DBSession.query(AnnotationFormField).join(AnnotationField).join(AnnotationItemFunctionalObject).join(ORM).all()


def __getAnnotationFields(type):
    if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
        from ..db_models.bioinformatics import Sequence as ORM
    elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
        from ..db_models.bioinformatics import MultipleSequenceAlignment as ORM
    elif type == 'phylotree' or type == 'phylotrees':
        from ..db_models.bioinformatics import PhylogeneticTree as ORM
    else:
        return None
    from ..db_models.sa_annotations import AnnotationField, AnnotationItemFunctionalObject
    _all = DBSession.query(AnnotationField).join(AnnotationItemFunctionalObject).join(ORM).all()
    form_fields = dict((i.form_field_id, i.form_field.name) for i in _all)

    import json
    from ..common import generate_json

    def _f(x):
        d = json.loads(generate_json(x))
        d['form_field'] = form_fields.get(x.form_field_id, '')
        return d

    return [_f(i) for i in _all]


def getFilterSchema(type, session):
    kwargs = {}
    kwargs['types'] = __getTypes(type)
    kwargs['cvterms'] = __getCvterms(type)
    kwargs['props'] = __getProps(type)
    kwargs['annotation_form_templates'] = __getAnnotationFormTemplates(type)
    kwargs['annotation_form_fields'] = __getAnnotationFormFields(type)
    kwargs['annotation_fields'] = __getAnnotationFields(type)
    kwargs['dbxref'] = __getDbxref(type)
    kwargs['organisms'] = __getOrganisms(type)
    kwargs['analyses'] = __getAnalyses(type)
    if kwargs.get('analyses'):
        kwargs['programs'] = [{'program': a.program} for a in kwargs.get('analyses')]
        kwargs['programs'] = [dict(t) for t in {tuple(p.items()) for p in kwargs['programs']}]
        kwargs['programversions'] = [{'program': a.program, 'programversion': a.programversion} for a in
                                     kwargs.get('analyses')]
        kwargs['programversions'] = [dict(t) for t in {tuple(p.items()) for p in kwargs['programversions']}]
        kwargs['algorithms'] = [{'algorithm': a.algorithm} for a in kwargs.get('analyses')]
        kwargs['algorithms'] = [dict(t) for t in {tuple(p.items()) for p in kwargs['algorithms']}]
        if type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
            kwargs['analyses'] = []
    elif type.startswith("annotation_form_"):
        from ..db_models.sa_annotations import AnnotationFormItem
        kwargs['standards'] = [{'label': i, 'value': i} for i, in
                               DBSession.query(AnnotationFormItem.standard).distinct(AnnotationFormItem.standard).all()]
    elif type == "geoprocesses_instances":  # Process instances
        kwargs["status"] = {
            'key': 'status',
            'type': 'customSelect',
            'templateOptions': {
                'nzMode': 'default',
                'label': 'Estado de la instancia:',
                'nzAllowClear': True,
                'nzShowSearch': True,
                'value': [],
                'options': [dict(label=f"{i[0]}", value=f"{i[1]}") for i in
                            [("Propuesta", "candidate"),
                             ("Planificada", "scheduled"),
                             ("Finalizada correctamente", "success"),
                             ("Finalizada con error", "error"),
                             ("Cancelada", "cancelled")]]
            }
        }
        case_studies = session.query(CaseStudy).all()
        kwargs["case_studies_"] = {
            'key': 'case_studies_',
            'type': 'customSelect',
            'templateOptions': {
                'nzMode': 'multiple',
                'label': 'Casos de estudio:',
                'nzAllowClear': True,
                'nzShowSearch': True,
                'value': [],
                'options': [dict(label=f"{i.name}", value=f"{i.id}") for i in case_studies]
            }
        }
        geoprocesses = session.query(CProcess).all()
        kwargs["for_geoprocesses"] = {
            'key': 'for_geoprocesses',
            'type': 'customSelect',
            'templateOptions': {
                'nzMode': 'multiple',
                'label': 'Instancia de geoprocesos:',
                'nzAllowClear': True,
                'nzShowSearch': True,
                'value': [],
                'options': [dict(label=f"{i.name}", value=f"{i.id}") for i in geoprocesses]
            }
        }
    elif type == "layers":  # Geographic Layers
        subjects = session.query(HierarchyNode).join(Hierarchy).\
            filter(Hierarchy.name == h_subjects_name).all()
        kwargs["subjects"] = {
            'key': 'subjects',
            'type': 'customSelect',
            'templateOptions': {
                'nzMode': 'multiple',
                'label': 'Temas:',
                'nzAllowClear': True,
                'nzShowSearch': True,
                'value': [],
                'options': [dict(label=f"{i.name}", value=f"{i.id}") for i in subjects]
            }
        }
        sources = session.query(HierarchyNode).join(Hierarchy).\
            filter(Hierarchy.name == h_sources_name).all()
        kwargs["sources"] = {
            'key': 'sources',
            'type': 'customSelect',
            'templateOptions': {
                'nzMode': 'default',
                'label': 'Fuente:',
                'nzAllowClear': True,
                'nzShowSearch': True,
                'value': [],
                'options': [dict(label=f"{i.name}", value=f"{i.id}") for i in sources]
            }
        }
        crs = session.query(HierarchyNode).join(Hierarchy).\
            filter(Hierarchy.name == h_crs_name).all()
        kwargs["crs"] = {
            'key': 'crs',
            'type': 'customSelect',
            'templateOptions': {
                'nzMode': 'default',
                'label': 'CRS:',
                'nzAllowClear': True,
                'nzShowSearch': True,
                'value': [],
                'options': [dict(label=f"{i.name}", value=f"{i.id}") for i in crs]
            }
        }
        case_studies = session.query(CaseStudy).all()
        kwargs["case_studies_"] = {
            'key': 'case_studies_',
            'type': 'customSelect',
            'templateOptions': {
                'nzMode': 'multiple',
                'label': 'Casos de estudio:',
                'nzAllowClear': True,
                'nzShowSearch': True,
                'value': [],
                'options': [dict(label=f"{i.name}", value=f"{i.id}") for i in case_studies]
            }
        }
        geoprocesses = session.query(CProcess).all()
        kwargs["used_as_input_of"] = {
            'key': 'used_as_input_of',
            'type': 'customSelect',
            'templateOptions': {
                'nzMode': 'multiple',
                'label': 'Usado en geoprocesos:',
                'nzAllowClear': True,
                'nzShowSearch': True,
                'value': [],
                'options': [dict(label=f"{i.name}", value=f"{i.id}") for i in geoprocesses]
            }
        }
        kwargs["resulting_from"] = {
            'key': 'resulting_from',
            'type': 'customSelect',
            'templateOptions': {
                'nzMode': 'multiple',
                'label': 'Resultante de geoprocesos:',
                'nzAllowClear': True,
                'nzShowSearch': True,
                'value': [],
                'options': [dict(label=f"{i.name}", value=f"{i.id}") for i in geoprocesses]
            }
        }

    from ..forms.filter_json import getJSONFilterSchema
    return getJSONFilterSchema(**kwargs)
