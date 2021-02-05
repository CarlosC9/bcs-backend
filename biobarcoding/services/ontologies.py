from biobarcoding.db_models import DBSessionChado as chado_session
from biobarcoding.rest import Issue, IType


def create_ontologies(name, definition = None, remote_url = None):
    issues = [Issue(IType.WARNING, 'CREATE ontologies: dummy completed')]
    content = { 'name':name, 'definition':definition, 'remote_url':remote_url }
    return issues, content, 200


def read_ontologies(ontology_id = None, name = None, definition = None, remote_url = None):
    content = {'ontology_id':ontology_id, 'name':name, 'definition':definition, 'remote_url':remote_url}
    content = {k:v for k,v in content.items() if v is not None}
    try:
        content = __get_query(cv_id=ontology_id, name=name, definition=definition)
        if ontology_id:
            content = content.first()
        else:
            content = content.all()
        issues, status = [Issue(IType.INFO, 'READ ontologies: The ontologies were successfully read.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, 'READ ontologies: The ontologies could not be read.')], 500
    return issues, content, status


def update_ontologies(ontology_id, name = None, definition = None, remote_url = None, input_file = None):
    issues = [Issue(IType.WARNING, 'UPDATE ontologies: dummy completed')]
    content = { 'name':name, 'definition':definition, 'remote_url':remote_url }
    return issues, content, 200


def delete_ontologies(ontology_id = None):
    issues = [Issue(IType.WARNING, 'DELETE ontologies: dummy completed')]
    content = { 'ontology_id':ontology_id }
    return issues, content, 200


def import_ontologies(input_file, format='obo', **kwargs):
    content = { 'input_file':input_file, 'format':format, **kwargs }
    content = {k:v for k,v in content.items() if v is not None}
    try:
        from flask import current_app
        cfg = current_app.config
        # f"""go2fmt.pl -p obo_text -w xml {input_file} | \
        #     go-apply-xslt oboxml_to_chadoxml - > {input_file}.xml"""
        # cmd = f"""go2chadoxml {input_file} > /tmp/{input_file}.chado.xml;
        #     stag-storenode.pl -d 'dbi:Pg:dbname={cfg['database']};host={cfg['host']};port=cfg['port']'
        #     --user {cfg['user']} --password {cfg['password']} /tmp/{input_file}.chado.xml"""
        import pronto
        import os
        namespace = pronto.Ontology(input_file).metadata.default_namespace
        if namespace:
            onto_name = f'-c {namespace}'
        else:
            onto_name = f'-c {os.path.basename(input_file)}'
        from biobarcoding.services import exec_cmds
        out, err = exec_cmds([
            f'''perl ./biobarcoding/services/perl_scripts/gmod_load_cvterms.pl\
                -H {cfg["CHADO_HOST"]}\
                -D {cfg["CHADO_DATABASE"]}\
                -r {cfg["CHADO_USER"]}\
                -p {cfg["CHADO_PASSWORD"]}\
                -d Pg -s null -u\
                {input_file}''',
            f'''perl ./biobarcoding/services/perl_scripts/gmod_make_cvtermpath.pl\
                -H {cfg["CHADO_HOST"]}\
                -D {cfg["CHADO_DATABASE"]}\
                -u {cfg["CHADO_USER"]}\
                -p {cfg["CHADO_PASSWORD"]}\
                -d Pg {onto_name}'''])
        issues, status = [Issue(IType.INFO, f'IMPORT ontologies: The {format} ontology were successfully imported.')], 200
    except Exception as e:
        print(e)
        issues, status = [Issue(IType.ERROR, f'IMPORT ontologies: The file {input_file} could not be imported.')], 500
    return issues, content, status


def export_ontologies(id=None, name=None, definition=None, remote_url=None, format='obo'):
    issues = [Issue(IType.WARNING, 'EXPORT ontologies: dummy completed')]
    content = { 'name':name, 'definition':definition, 'remote_url':remote_url, 'format':format }
    return issues, content, 200


def __get_query(**kwargs):
    from biobarcoding.db_models.chado import Cv
    result = chado_session.query(Cv)
    if kwargs.get('ontology_id'):
        result = result.filter(Cv.cv_id==kwargs.get('ontology_id'))
    if kwargs.get('name'):
        result = result.filter(Cv.name==kwargs.get('name'))
    if kwargs.get('definition'):
        result = result.filter(Cv.definition==kwargs.get('definition'))
    return result


# TODO:
#  featureprop, analysisprop, phylotreeprop ?
#  phylotree ?
def read_cvterms(cv_id=None, cvterm_id=None, feature_id=None, analysis_id=None, phylotree_id=None):
    content = {'cv_id':cv_id, 'cvterm_id':cvterm_id, 'feature_id':feature_id, 'analysis_id':analysis_id, 'phylotree_id':phylotree_id}
    content = {k:v for k,v in content.items() if v is not None}
    try:
        from biobarcoding.db_models.chado import Cvterm
        result = chado_session.query(Cvterm)
        if cvterm_id:
            result = result.filter(Cvterm.cvterm_id==cvterm_id)
        if cv_id:
            result = result.filter(Cvterm.cv_id==cv_id)
        if feature_id:
            from biobarcoding.db_models.chado import FeatureCvterm
            cv_ids = chado_session.query(FeatureCvterm.cvterm_id)\
                .filter(FeatureCvterm.feature_id==feature_id)
            result = result.filter(Cvterm.cv_id.in_(cv_ids))
        if analysis_id:
            from biobarcoding.db_models.chado import AnalysisCvterm
            cv_ids = chado_session.query(AnalysisCvterm.cvterm_id)\
                .filter(AnalysisCvterm.analysis_id==analysis_id)
            result = result.filter(Cvterm.cv_id.in_(cv_ids))
        if phylotree_id:
            # from biobarcoding.db_models.chado import Phylotree
            # cv_id = chado_session.query(Phylotree.type_id)\
            #     .filter(Phylotree.phylotree_id==phylotree_id)
            # result = result.filter(Cvterm.cv_id==cv_id)
            pass
        if cvterm_id:
            content = result.first()
        else:
            content = result.all()
        issues, status = [Issue(IType.ERROR, f'READ ontology terms: The cvterms were successfully read.')], 200
    except Exception as e:
        issues, status = [Issue(IType.ERROR, f'READ ontology terms: The cvterms could not be read.')], 500
    return issues, content, status
