from flask import Blueprint, request
from flask.views import MethodView

from biobarcoding.rest.bos import BioObjAPI

bp_metadata = Blueprint('bp_metadata', __name__)

from biobarcoding.authentication import bcs_session
from biobarcoding.rest import bcs_api_base, ResponseObject

metadata = ['taxonomies', 'organisms', 'ontologies', 'analyses' ]

class MetadataAPI(BioObjAPI):
    """
    Taxonomies Resource
    """

    @bcs_session(read_only=True)
    def get(self, id=None, format=None):
        return super(MetadataAPI, self).get(self.getMetadata(), id, format)

    @bcs_session()
    def post(self, format=None):
        return super(MetadataAPI, self).post(self.getMetadata(), format)

    @bcs_session()
    def put(self, id, format=None):
        return super(MetadataAPI, self).put(self.getMetadata(), id, format)

    @bcs_session()
    def delete(self, id=None, format=None):
        return super(MetadataAPI, self).delete(self.getMetadata(), id, format)

    def getMetadata(self):
        for meta in metadata:
            # if request.url_rule.rule.split(bcs_api_base+'/')[1].startswith(meta):
            if meta in request.url_rule.rule:
                return meta


metadata_view = MetadataAPI.as_view('api_metadata')

for meta in metadata:
    bp_metadata.add_url_rule(
        bcs_api_base + f'/{meta}/',
        view_func=metadata_view,
        methods=['GET', 'POST', 'DELETE']
    )
    bp_metadata.add_url_rule(
        bcs_api_base + f'/{meta}/<string:id>',
        view_func=metadata_view,
        methods=['GET', 'PUT', 'DELETE']
    )
    bp_metadata.add_url_rule(
        bcs_api_base + f'/{meta}.<string:format>',
        view_func=metadata_view,
        methods=['GET', 'POST']
    )
    bp_metadata.add_url_rule(
        bcs_api_base + f'/{meta}/<string:id>.<string:format>',
        view_func=metadata_view,
        methods=['GET', 'PUT', 'DELETE']
    )


class CvtermsAPI(BioObjAPI):
    """
    Cvterms Resource
    """

    @bcs_session(read_only=True)
    def get(self, cv_id=None, cvterm_id=None):
        print(f'GET {request.path}\nGetting ontology terms {id}')
        self._check_data()
        from biobarcoding.services.ontologies import read_cvterms
        issues, content, status = read_cvterms(cv_id, cvterm_id, **self.kwargs)
        return ResponseObject(content=content, issues=issues, status=status).get_response()


cvterms_view = CvtermsAPI.as_view('api_cvterms')
bp_metadata.add_url_rule(
    bcs_api_base + '/ontologies/terms/',
    view_func=cvterms_view,
    methods=['GET']
)
bp_metadata.add_url_rule(
    bcs_api_base + '/ontologies/<int:cv_id>/terms/',
    view_func=cvterms_view,
    methods=['GET']
)
bp_metadata.add_url_rule(
    bcs_api_base + '/ontologies/terms/<int:cvterm_id>',
    view_func=cvterms_view,
    methods=['GET']
)
bp_metadata.add_url_rule(
    bcs_api_base + '/ontologies/<int:cv_id>/terms/<int:cvterm_id>',
    view_func=cvterms_view,
    methods=['GET']
)
