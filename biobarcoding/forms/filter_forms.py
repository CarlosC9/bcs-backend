from biobarcoding.db_models import DBSessionChado as chado_session


def __getTypes(type):
  if type == 'sequence' or type == 'sequences':
    from biobarcoding.db_models.chado import Feature as ORM
  elif type == 'phylotree' or type == 'phylotrees':
    from biobarcoding.db_models.chado import Phylotree as ORM
  else:
    return None
  ids = chado_session.query(ORM.type_id).distinct(ORM.type_id).all()
  from biobarcoding.services.ontologies import __get_cvterm as read_cvterms
  ids = read_cvterms(filter=[{'cvterm_id': {'op': 'in', 'unary': ids}}]).all()
  return ids


def __getCvterms(type):
  if type == 'sequence' or type == 'sequences':
    from biobarcoding.db_models.chado import FeatureCvterm as ORM
  elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
    from biobarcoding.db_models.chado import AnalysisCvterm as ORM
  else:
    return None
  ids = chado_session.query(ORM.cvterm_id).distinct(ORM.cvterm_id).all()
  from biobarcoding.services.ontologies import __get_cvterm as read_cvterms
  ids = read_cvterms(filter=[{'cvterm_id': {'op': 'in', 'unary': ids}}]).all()
  return ids


def __getProps(type):
  if type == 'sequence' or type == 'sequences':
    from biobarcoding.db_models.chado import Featureprop as ORM
  elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
    from biobarcoding.db_models.chado import Analysisprop as ORM
  elif type == 'phylotree' or type == 'phylotrees':
    from biobarcoding.db_models.chado import Phylotreeprop as ORM
  else:
    return None
  ids = chado_session.query(ORM.type_id).distinct(ORM.type_id).all()
  from biobarcoding.services.ontologies import __get_cvterm as read_cvterms
  ids = read_cvterms(filter=[{'cvterm_id': {'op': 'in', 'unary': ids}}]).all()
  return ids


def __getDbxref(type):
  if type == 'sequence' or type == 'sequences':
    from biobarcoding.db_models.chado import Feature as ORM
  elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
    from biobarcoding.db_models.chado import AnalysisDbxref as ORM
  elif type == 'phylotree' or type == 'phylotrees':
    from biobarcoding.db_models.chado import Phylotree as ORM
  else:
    return None
  ids = chado_session.query(ORM.dbxref_id).distinct(ORM.dbxref_id).all()
  from biobarcoding.db_models.chado import Dbxref
  ids = chado_session.query(Dbxref).filter(Dbxref.dbxref_id.in_(ids)).all()
  return ids


def __getOrganisms(type):
  from biobarcoding.db_models.chado import Feature as ORM
  if type == 'phylotree' or type == 'phylotrees':
    # Phylonode > Feature
    from biobarcoding.db_models.chado import Phylonode
    ids = chado_session.query(Phylonode.feature_id).distinct(Phylonode.feature_id).all()
    ids = chado_session.query(ORM.organism_id).filter(ORM.feature_id.in_(ids))
  # if type == 'taxonomy' or type == 'taxonomies':
    # Phylonode > PhylonodeOrganism
    # from biobarcoding.db_models.chado import PhylonodeOrganism
    # ids += chado_session.query(PhylonodeOrganism.organism_id).distinct(PhylonodeOrganism.organism_id).all()
  elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
    # AnalysisFeature > Feature
    from biobarcoding.db_models.chado import AnalysisFeature
    ids = chado_session.query(AnalysisFeature.feature_id).distinct(AnalysisFeature.feature_id).all()
    ids = chado_session.query(ORM.organism_id).filter(ORM.feature_id.in_(ids))
  elif type == 'sequence' or type == 'sequences':
    ids = chado_session.query(ORM.organism_id)
  else:
    return None
  from biobarcoding.db_models.chado import Organism
  ids = ids.distinct(ORM.organism_id).all()
  ids = chado_session.query(Organism).filter(Organism.organism_id.in_(ids)).all()
  return ids


def __getAnalyses(type):
  if type == 'sequence' or type == 'sequences':
    from biobarcoding.db_models.chado import AnalysisFeature as ORM
  elif type == 'phylotree' or type == 'phylotrees':
    from biobarcoding.db_models.chado import Phylotree as ORM
  else:
    return None
  ids = chado_session.query(ORM.analysis_id).distinct(ORM.analysis_id).all()
  from biobarcoding.db_models.chado import Analysis
  ids = chado_session.query(Analysis).filter(Analysis.analysis_id.in_(ids)).all()
  return ids


def getFilterSchema(type):
  kwargs = {}
  kwargs['types'] = __getTypes(type)
  kwargs['cvterms'] = __getCvterms(type)
  kwargs['props'] = __getProps(type)
  kwargs['dbxref'] = __getDbxref(type)
  kwargs['organisms'] = __getOrganisms(type)
  kwargs['analyses'] = __getAnalyses(type)
  from biobarcoding.forms.filter_json import getJSONFilterSchema
  return getJSONFilterSchema(**kwargs)
