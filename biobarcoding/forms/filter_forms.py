from biobarcoding.db_models import DBSessionChado as chado_session


def __getTypes(type):
  if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
    from biobarcoding.db_models.chado import Feature as ORM
  elif type == 'phylotree' or type == 'phylotrees':
    from biobarcoding.db_models.chado import Phylotree as ORM
  elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
    from biobarcoding.db_models.chado import Stock as ORM
  elif type == 'collection' or type == 'collections' or type == 'stockcollection' or type == 'stockcollections':
    from biobarcoding.db_models.chado import Stockcollection as ORM
  else:
    return None
  ids = chado_session.query(ORM.type_id).distinct(ORM.type_id).all()
  from biobarcoding.services.ontologies import __get_cvterm as read_cvterms
  ids = read_cvterms(filter=[{'cvterm_id': {'op': 'in', 'unary': ids}}]).all()
  return ids


def __getCvterms(type):
  if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
    from biobarcoding.db_models.chado import FeatureCvterm as ORM
  elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
    from biobarcoding.db_models.chado import AnalysisCvterm as ORM
  elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
    from biobarcoding.db_models.chado import StockCvterm as ORM
  else:
    return None
  ids = chado_session.query(ORM.cvterm_id).distinct(ORM.cvterm_id).all()
  from biobarcoding.services.ontologies import __get_cvterm as read_cvterms
  ids = read_cvterms(filter=[{'cvterm_id': {'op': 'in', 'unary': ids}}]).all()
  return ids


def __getProps(type):
  if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
    from biobarcoding.db_models.chado import Featureprop as ORM
  elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
    from biobarcoding.db_models.chado import Analysisprop as ORM
  elif type == 'phylotree' or type == 'phylotrees':
    from biobarcoding.db_models.chado import Phylotreeprop as ORM
  elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
    from biobarcoding.db_models.chado import Stockprop as ORM
  elif type == 'collection' or type == 'collections' or type == 'stockcollection' or type == 'stockcollections':
    from biobarcoding.db_models.chado import Stockcollectionprop as ORM
  else:
    return None
  ids = chado_session.query(ORM.type_id).distinct(ORM.type_id).all()
  from biobarcoding.services.ontologies import __get_cvterm as read_cvterms
  ids = read_cvterms(filter=[{'cvterm_id': {'op': 'in', 'unary': ids}}]).all()
  return ids


def __getDbxref(type):
  if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
    from biobarcoding.db_models.chado import Feature as ORM
  elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
    from biobarcoding.db_models.chado import AnalysisDbxref as ORM
  elif type == 'phylotree' or type == 'phylotrees':
    from biobarcoding.db_models.chado import Phylotree as ORM
  elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
    from biobarcoding.db_models.chado import Stock as ORM
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
  elif type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
    ids = chado_session.query(ORM.organism_id)

  elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
    from biobarcoding.db_models.chado import Stock as ORM
    ids = chado_session.query(ORM.organism_id)
  elif type == 'collection' or type == 'collections' or type == 'stockcollection' or type == 'stockcollections':
    from biobarcoding.db_models.chado import Stock as ORM, StockcollectionStock
    ids = chado_session.query(StockcollectionStock.stock_id).distinct(StockcollectionStock.stock_id).all()
    ids = chado_session.query(ORM.organism_id).filter(ORM.stock_id.in_(ids))
  else:
    return None
  from biobarcoding.db_models.chado import Organism
  ids = ids.distinct(ORM.organism_id).all()
  ids = chado_session.query(Organism).filter(Organism.organism_id.in_(ids)).all()
  return ids


def __getAnalyses(type):
  if type == 'sequence' or type == 'sequences' or type == 'feature' or type == 'features':
    from biobarcoding.db_models.chado import AnalysisFeature as ORM
  elif type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
    from biobarcoding.db_models.chado import Analysis as ORM
  elif type == 'phylotree' or type == 'phylotrees':
    from biobarcoding.db_models.chado import Phylotree as ORM
  # elif type == 'individual' or type == 'individuals' or type == 'stock' or type == 'stocks':
  #   from biobarcoding.db_models.chado import StockFeature
  #   ids = chado_session.query(StockFeature.feature_id).distinct(StockFeature.feature_id)
  # elif type == 'collection' or type == 'collections' or type == 'stockcollection' or type == 'stockcollections':
  #   from biobarcoding.db_models.chado import Stock as ORM, StockcollectionStock
  #   ids = chado_session.query(StockcollectionStock.stock_id).distinct(StockcollectionStock.stock_id).all()
  #   ids = chado_session.query(ORM.organism_id).filter(ORM.stock_id.in_(ids))
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
  if kwargs.get('analyses'):
    kwargs['programs'] = [ {'program':a.program} for a in kwargs.get('analyses') ]
    kwargs['programs'] = [dict(t) for t in {tuple(p.items()) for p in kwargs['programs']}]
    kwargs['programversions'] = [ {'program':a.program, 'programversion':a.programversion} for a in kwargs.get('analyses') ]
    kwargs['programversions'] = [dict(t) for t in {tuple(p.items()) for p in kwargs['programversions']}]
    kwargs['algorithms'] = [ {'algorithm':a.algorithm} for a in kwargs.get('analyses') ]
    kwargs['algorithms'] = [dict(t) for t in {tuple(p.items()) for p in kwargs['algorithms']}]
    if type == 'alignment' or type == 'analysis' or type == 'alignments' or type == 'analyses':
      kwargs['analyses']=[]
  from biobarcoding.forms.filter_json import getJSONFilterSchema
  return getJSONFilterSchema(**kwargs)
