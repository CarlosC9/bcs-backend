
datetime = [{
    'key': 'added-from',
    'type': 'input',
    'templateOptions': {
        'type': 'date',
        'label': 'Añadido después de',
    },
}, {
    'key': 'added-to',
    'type': 'input',
    'templateOptions': {
        'type': 'date',
        'label': 'Añadido antes de',
    }
}, {
    'key': 'lastmodified-from',
    'type': 'input',
    'templateOptions': {
        'type': 'date',
        'label': 'Modificado después de',
    },
}, {
    'key': 'lastmodified-to',
    'type': 'input',
    'templateOptions': {
        'type': 'date',
        'label': 'Modificado antes de',
    }
}]


def getJSONFilterSchema(**kwargs):
    schema = []
    for key in kwargs:
        # if there are values for the filter and the filter itself
        if kwargs[key]:
            schema.append(__getJSONFilter(key, kwargs[key]))
    return [ x for x in schema if x ] + datetime


formly = {'types':
               {'key': 'type_id',
                'type': 'select',
                'templateOptions': {
                    'label': 'Tipo',
                    'placeholder': 'Tipo de elementos',
                    'multiple': True,
                    'options': [],
                    'valueProp': 'cvterm_id',
                    'labelProp': 'name',
                }},
           'cvterms':
               {'key': 'cvterm_id',
                'type': 'select',
                'templateOptions': {
                    'label': 'Términos',
                    'placeholder': 'Términos asociados',
                    'multiple': True,
                    'options': [],
                    'valueProp': 'cvterm_id',
                    'labelProp': 'name',
                }},
           'props':
           # TODO: prop-value-selector. how?
               {'key': 'prop_cvterm_id',
                'type': 'select',
                'templateOptions': {
                    'label': 'Propiedades',
                    'placeholder': 'Propiedades asociadas',
                    'multiple': True,
                    'options': [],
                    'valueProp': 'cvterm_id',
                    'labelProp': 'name',
                }},
           'dbxrefs':
               {'key': 'dbxref_id',
                'type': 'select',
                'templateOptions': {
                    'label': 'Referencias',
                    'placeholder': 'Referencias externas',
                    'multiple': True,
                    'options': [],
                    'valueProp': 'dbxref_id',
                    'labelProp': 'accession',
                }},
           'organisms':
           # TODO: avoid clearing search after choosing
               {'key': 'organism_id',
                'type': 'select',
                'templateOptions': {
                    'label': 'Taxones',
                    'placeholder': 'Taxones de interés',
                    'multiple': True,
                    'options': [],
                    'groupProp': 'genus',
                    'valueProp': 'organism_id',
                    'labelProp': 'species',
                }},
           'analyses':
               {'key': 'analysis_id',
                'type': 'select',
                'templateOptions': {
                    'label': 'Análisis',
                    'placeholder': 'Análisis de interés',
                    'multiple': True,
                    'options': [],
                    'groupProp': 'program',
                    'valueProp': 'analysis_id',
                    'labelProp': 'name',
                }},
           'programs':
               {'key': 'program',
                'type': 'select',
                'templateOptions': {
                    'label': 'Programa',
                    'placeholder': 'Programa de interés',
                    'multiple': True,
                    'options': [],
                    'valueProp': 'program',
                    'labelProp': 'program',
                }},
           'programversions':
               {'key': 'programversion',
                'type': 'select',
                'templateOptions': {
                    'label': 'Programa versión',
                    'placeholder': 'Versiones de interés',
                    'multiple': True,
                    'options': [],
                    'groupProp': 'program',
                    'valueProp': 'programversion',
                    'labelProp': 'programversion',
                }},
           'algorithms':
               {'key': 'algorithm',
                'type': 'select',
                'templateOptions': {
                    'label': 'Algoritmo',
                    'placeholder': 'Algoritmo de interés',
                    'multiple': True,
                    'options': [],
                    'valueProp': 'algorithm',
                    'labelProp': 'algorithm',
                }}
          }


def __getJSONFilter(type, value):
    if not type in formly:
        return None
    filter = formly[type].copy()
    filter['templateOptions']['options'] = value
    return filter
