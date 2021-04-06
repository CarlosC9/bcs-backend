# params
# [
#   { label: '', value: '' },
#   { label: '', value: '' }
# ]
def getJSONFilterSchema(type, types=[], cvterms=[], props=[], dbxref=[], organisms=[], analyses=[]):
    schemas = {
      'default': [ {
        'key': 'types',
        'type': 'select',
        'templateOptions': {
          'label': 'Tipo',
          'placeholder': 'Tipo de elementos',
          'multiple': True,
          'options': types,
        },
      }, {
        'key': 'cvterms',
        'type': 'select',
        'templateOptions': {
          'label': 'Términos',
          'placeholder': 'Términos asociados',
          'multiple': True,
          'options': cvterms,
        },
      }, {
        'key': 'properties',
        'type': 'select',
        'templateOptions': {
          'label': 'Propiedades',
          'placeholder': 'Propiedades asociadas',
          'multiple': True,
          'options': props,
        },
        # TODO: prop-value-selector. how?
      },{
        'key': 'added-from',
        'type': 'input',
        'templateOptions': {
          'type': 'date',
          'label': 'Añadido antes de',
        },
      },{
        'key': 'added-to',
        'type': 'input',
        'templateOptions': {
          'type': 'date',
          'label': 'Añadido después de',
        }
      },{
        'key': 'lastmodified-from',
        'type': 'input',
        'templateOptions': {
          'type': 'date',
          'label': 'Modificado antes de',
        },
      },{
        'key': 'lastmodified-to',
        'type': 'input',
        'templateOptions': {
          'type': 'date',
          'label': 'Modificado después de',
        }
      },{
        'key': 'dbxref',
        'type': 'select',
        'templateOptions': {
          'label': 'Referencias',
          'placeholder': 'Referencias externas',
          'multiple': True,
          'options': dbxref,
        }
      } ],

      'sequence': [ {
        # TODO: avoid clearing search after choosing
        'key': 'organisms',
        'type': 'select',
        'templateOptions': {
          'label': 'Taxones',
          'placeholder': 'Taxones de interés',
          'multiple': True,
          'options': organisms,
          'groupProp': 'genus',
          'valueProp': 'organism_id',
          'labelProp': 'species',

        },
      },{
        'key': 'analyses',
        'type': 'select',
        'templateOptions': {
          'label': 'Análisis',
          'placeholder': 'Análisis de interés',
          'multiple': True,
          'options': analyses,
          'groupProp': 'program',
          'valueProp': 'analysis_id',
          'labelProp': 'name',
        },
      } ],

      'alignment': [],

      'phylotree': [],

      'analysis': [],
      'taxonomy': [],
      'organism': [],
      'ontology': [],
      'cvterm': [],
    }
    return schemas['default'] + schemas[type]