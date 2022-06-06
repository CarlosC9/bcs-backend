from ..db_models.chado import Cvterm
from ..db_models.core import CaseStudy, CProcess
from ..db_models.hierarchies import HierarchyNode, Hierarchy
from ..rest import h_subjects_name, h_sources_name, h_crs_name
from ..db_models import DBSession, DBSessionChado


def __getTypes(subject):
    if subject == 'sequence' or subject == 'sequences' or subject == 'feature' or subject == 'features':
        from ..db_models.chado import Feature as ORM
    elif subject == 'phylotree' or subject == 'phylotrees':
        from ..db_models.chado import Phylotree as ORM
    elif subject == 'individual' or subject == 'individuals' or subject == 'stock' or subject == 'stocks':
        from ..db_models.chado import Stock as ORM
    else:
        return None
    ids = DBSessionChado.query(ORM.type_id).distinct(ORM.type_id).subquery()
    return DBSessionChado.query(Cvterm).filter(Cvterm.cvterm_id.in_(ids)).all()


def __getCvterms(subject):
    if subject == 'sequence' or subject == 'sequences' or subject == 'feature' or subject == 'features':
        from ..db_models.chado import FeatureCvterm as ORM
    elif subject == 'alignment' or subject == 'analysis' or subject == 'alignments' or subject == 'analyses':
        from ..db_models.chado import AnalysisCvterm as ORM
    elif subject == 'individual' or subject == 'individuals' or subject == 'stock' or subject == 'stocks':
        from ..db_models.chado import StockCvterm as ORM
    elif subject == 'annotation_form_template' or subject == 'annotation_form_templates' \
            or subject == 'annotation_form_field' or subject == 'annotation_form_fields':
        from ..db_models.sa_annotations import AnnotationFormItem as ORM
        ids = DBSession.query(ORM.cvterm_id).distinct(ORM.cvterm_id).all()
        return DBSessionChado.query(Cvterm).filter(Cvterm.cvterm_id.in_(ids)).all()
    else:
        return None
    ids = DBSessionChado.query(ORM.cvterm_id).distinct(ORM.cvterm_id).subquery()
    return DBSessionChado.query(Cvterm).filter(Cvterm.cvterm_id.in_(ids)).all()


def __getProps(subject):
    if subject == 'sequence' or subject == 'sequences' or subject == 'feature' or subject == 'features':
        from ..db_models.chado import Featureprop as ORM
    elif subject == 'alignment' or subject == 'analysis' or subject == 'alignments' or subject == 'analyses':
        from ..db_models.chado import Analysisprop as ORM
    elif subject == 'phylotree' or subject == 'phylotrees':
        from ..db_models.chado import Phylotreeprop as ORM
    elif subject == 'individual' or subject == 'individuals' or subject == 'stock' or subject == 'stocks':
        from ..db_models.chado import Stockprop as ORM
    else:
        return None
    ids = DBSessionChado.query(ORM.type_id).distinct(ORM.type_id).subquery()
    return DBSessionChado.query(Cvterm).filter(Cvterm.cvterm_id.in_(ids)).all()


def __getDbxref(subject):
    from ..db_models.chado import Dbxref
    if subject == 'sequence' or subject == 'sequences' or subject == 'feature' or subject == 'features':
        from ..db_models.chado import Feature as ORM
    elif subject == 'alignment' or subject == 'analysis' or subject == 'alignments' or subject == 'analyses':
        from ..db_models.chado import AnalysisDbxref as ORM
    elif subject == 'phylotree' or subject == 'phylotrees':
        from ..db_models.chado import Phylotree as ORM
    elif subject == 'individual' or subject == 'individuals' or subject == 'stock' or subject == 'stocks':
        from ..db_models.chado import Stock as ORM
    elif subject == 'annotation_form_template' or subject == 'annotation_form_templates' \
            or subject == 'annotation_form_field' or subject == 'annotation_form_fields':
        from ..db_models.sa_annotations import AnnotationFormItem as ORM
        ids = DBSession.query(ORM.dbxref_id).distinct(ORM.dbxref_id).all()
        return DBSessionChado.query(Dbxref).filter(Dbxref.dbxref_id.in_(ids)).all()
    else:
        return None
    ids = DBSessionChado.query(ORM.dbxref_id).distinct(ORM.dbxref_id).subquery()
    return DBSessionChado.query(Dbxref).filter(Dbxref.dbxref_id.in_(ids)).all()


def __getOrganisms(subject):
    from ..db_models.chado import Feature as ORM
    if subject == 'sequence' or subject == 'sequences' or subject == 'feature' or subject == 'features':
        ids = DBSessionChado.query(ORM.organism_id)
    elif subject == 'alignment' or subject == 'analysis' or subject == 'alignments' or subject == 'analyses':
        # AnalysisFeature > Feature
        from ..db_models.chado import AnalysisFeature
        ids = DBSessionChado.query(AnalysisFeature.feature_id).subquery()
        ids = DBSessionChado.query(ORM.organism_id).filter(ORM.feature_id.in_(ids))
    elif subject == 'phylotree' or subject == 'phylotrees':
        # Phylonode > Feature
        from ..db_models.chado import Phylonode
        ids = DBSessionChado.query(Phylonode.feature_id).subquery()
        ids = DBSessionChado.query(ORM.organism_id).filter(ORM.feature_id.in_(ids))
    elif subject == 'individual' or subject == 'individuals' or subject == 'stock' or subject == 'stocks':
        from ..db_models.chado import Stock as ORM
        ids = DBSessionChado.query(ORM.organism_id)
    elif subject == 'taxonomy' or subject == 'taxonomies':
        from ..db_models.chado import PhylonodeOrganism as ORM
        ids = DBSessionChado.query(ORM.organism_id)
    else:
        return None
    ids = ids.distinct(ORM.organism_id).all()
    from ..services.bio.meta.organisms import Service as OrgService
    return OrgService().read(filter={'organism_id': ids})[0]    # info attachment needed


def __getAnalyses(subject):
    if subject == 'sequence' or subject == 'sequences' or subject == 'feature' or subject == 'features':
        from ..db_models.chado import AnalysisFeature as ORM
    elif subject == 'alignment' or subject == 'analysis' or subject == 'alignments' or subject == 'analyses':
        from ..db_models.chado import Analysis as ORM
    elif subject == 'phylotree' or subject == 'phylotrees':
        from ..db_models.chado import Phylotree as ORM
    else:
        return None
    ids = DBSessionChado.query(ORM.analysis_id).distinct(ORM.analysis_id).subquery()
    from ..db_models.chado import Analysis
    return DBSessionChado.query(Analysis).filter(Analysis.analysis_id.in_(ids)).all()


def __getAnnotationFormTemplates(subject):
    if subject == 'sequence' or subject == 'sequences' or subject == 'feature' or subject == 'features':
        from ..db_models.bioinformatics import Sequence as ORM
    elif subject == 'alignment' or subject == 'analysis' or subject == 'alignments' or subject == 'analyses':
        from ..db_models.bioinformatics import MultipleSequenceAlignment as ORM
    elif subject == 'phylotree' or subject == 'phylotrees':
        from ..db_models.bioinformatics import PhylogeneticTree as ORM
    elif subject == 'annotation_form_field' or subject == 'annotation_form_fields':
        from ..db_models.sa_annotations import AnnotationFormTemplate, AnnotationFormTemplateField as ORM
        return DBSession.query(AnnotationFormTemplate).join(ORM).all()
    else:
        return None
    from ..db_models.sa_annotations import AnnotationFormTemplate, AnnotationTemplate, AnnotationItemFunctionalObject
    return DBSession.query(AnnotationFormTemplate).join(AnnotationTemplate).join(AnnotationItemFunctionalObject).join(ORM).all()


def __getAnnotationFormFields(subject):
    if subject == 'sequence' or subject == 'sequences' or subject == 'feature' or subject == 'features':
        from ..db_models.bioinformatics import Sequence as ORM
    elif subject == 'alignment' or subject == 'analysis' or subject == 'alignments' or subject == 'analyses':
        from ..db_models.bioinformatics import MultipleSequenceAlignment as ORM
    elif subject == 'phylotree' or subject == 'phylotrees':
        from ..db_models.bioinformatics import PhylogeneticTree as ORM
    elif subject == 'annotation_form_template' or subject == 'annotation_form_templates':
        from ..db_models.sa_annotations import AnnotationFormField, AnnotationFormTemplateField as ORM
        return DBSession.query(AnnotationFormField).join(ORM).all()
    else:
        return None
    from ..db_models.sa_annotations import AnnotationFormField, AnnotationField, AnnotationItemFunctionalObject
    return DBSession.query(AnnotationFormField).join(AnnotationField).join(AnnotationItemFunctionalObject).join(ORM).all()


def __getAnnotationFields(subject):
    if subject == 'sequence' or subject == 'sequences' or subject == 'feature' or subject == 'features':
        from ..db_models.bioinformatics import Sequence as ORM
    elif subject == 'alignment' or subject == 'analysis' or subject == 'alignments' or subject == 'analyses':
        from ..db_models.bioinformatics import MultipleSequenceAlignment as ORM
    elif subject == 'phylotree' or subject == 'phylotrees':
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


def getFilterSchema(subject, session):
    kwargs = {}
    kwargs['types'] = __getTypes(subject)
    kwargs['cvterms'] = __getCvterms(subject)
    kwargs['props'] = __getProps(subject)
    kwargs['annotation_form_templates'] = __getAnnotationFormTemplates(subject)
    kwargs['annotation_form_fields'] = __getAnnotationFormFields(subject)
    kwargs['annotation_fields'] = __getAnnotationFields(subject)
    kwargs['dbxref'] = __getDbxref(subject)
    kwargs['organisms'] = __getOrganisms(subject)
    kwargs['analyses'] = __getAnalyses(subject)
    if kwargs.get('analyses'):
        kwargs['programs'] = [{'program': a.program} for a in kwargs.get('analyses')]
        kwargs['programs'] = [dict(t) for t in {tuple(p.items()) for p in kwargs['programs']}]
        kwargs['programversions'] = [{'program': a.program, 'programversion': a.programversion} for a in
                                     kwargs.get('analyses')]
        kwargs['programversions'] = [dict(t) for t in {tuple(p.items()) for p in kwargs['programversions']}]
        kwargs['algorithms'] = [{'algorithm': a.algorithm} for a in kwargs.get('analyses')]
        kwargs['algorithms'] = [dict(t) for t in {tuple(p.items()) for p in kwargs['algorithms']}]
        if subject == 'alignment' or subject == 'analysis' or subject == 'alignments' or subject == 'analyses':
            kwargs['analyses'] = []
    elif subject.startswith("annotation_form_"):
        from ..db_models.sa_annotations import AnnotationFormItem
        kwargs['standards'] = [{'label': i, 'value': i} for i, in
                               DBSession.query(AnnotationFormItem.standard).distinct(AnnotationFormItem.standard).all()]
    elif subject == "geoprocesses_instances":  # Process instances
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
    elif subject == "layers":  # Geographic Layers
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
    return getJSONFilterSchema(subject, **kwargs)
