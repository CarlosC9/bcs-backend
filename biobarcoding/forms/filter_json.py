from ..services.bio.meta.ontologies import get_type_id

datetime_filter_fields = [{
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
          'standards':
              {'key': 'standard',
               'type': 'select',
               'templateOptions': {
                   'label': 'Estándares',
                   'placeholder': 'Estándares asociados',
                   'multiple': True,
                   'options': [],
               }},
          'annotation_form_templates':
              {'key': 'annotation_form_template_id',
               'type': 'select',
               'templateOptions': {
                   'label': 'Plantillas',
                   'placeholder': 'Plantillas asociados',
                   'multiple': True,
                   'options': [],
                   # 'groupProp': 'type',
                   'groupProp': 'standard',
                   'valueProp': 'id',
                   'labelProp': 'name',
               }},
          'annotation_form_fields':
              {'key': 'annotation_form_field_id',
               'type': 'select',
               'templateOptions': {
                   'label': 'Campos',
                   'placeholder': 'Campos asociados',
                   'multiple': True,
                   'options': [],
                   # 'groupProp': 'type',
                   'groupProp': 'standard',
                   'valueProp': 'id',
                   'labelProp': 'name',
               }},
          'annotation_fields':
              {'key': 'annotation_field_id',
               'type': 'select',
               'templateOptions': {
                   'label': 'Propiedades',
                   'placeholder': 'Propiedades asociadas',
                   'multiple': True,
                   'options': [],
                   'groupProp': 'form_field',
                   'valueProp': 'id',
                   'labelProp': 'value',
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
                   'labelProp': 'canonical_name',
               }},
          'genus':
              {'key': 'genus',
               'type': 'select',
               'templateOptions': {
                   'label': 'Género',
                   'placeholder': 'Géneros de interés',
                   'multiple': True,
                   'options': [],
                   'valueProp': 'genus',
                   'labelProp': 'genus',
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


def getJSONFilterSchema(subject=None, **kwargs):
    schema = []
    for key in kwargs:
        # if there are values for the filter and the filter itself
        if key not in formly:
            schema.append(kwargs[key])
        elif kwargs.get(key):
            _ = __getJSONFilter(key, kwargs[key])
            if subject.startswith('sequence') and key == 'types':
                _['defaultValue'] = [get_type_id('sequence')]
            schema.append(_)
    return [x for x in schema if x] + datetime_filter_fields

