def __getTypes(type):
  # import service, get objs, and filter those that are bind to some
  pass


def __getCvterms(type):
  # import service, get objs, and filter those that are bind to some
  pass


def __getProps(type):
  # import service, get objs, and filter those that are bind to some
  pass


def __getDbxref(type):
  # import service, get objs, and filter those that are bind to some
  pass


def __getOrganisms(type):
  # import service, get objs, and filter those that are bind to some
  pass


def __getAnalyses(type):
  # import service, get objs, and filter those that are bind to some
  pass


def getFilterSchema(type):
  # types = __getTypes(type)
  # cvterms = __getCvterms(type)
  # props = __getProps(type)
  # dbxref = __getDbxref(type)
  # organisms = __getOrganisms(type)
  # analyses = __getAnalyses(type)
  # return getFilterSchema(type, types, cvterms, props, dbxref, organisms, analyses)
  from biobarcoding.forms.filter_json import getJSONFilterSchema
  return getJSONFilterSchema(type)
