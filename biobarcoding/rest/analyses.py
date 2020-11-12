from flask import Blueprint

bp_analyses = Blueprint('bp_analyses', __name__)

from flask import request, make_response, jsonify
from flask.views import MethodView

from biobarcoding.rest import bcs_api_base


class AnalysesAPI(MethodView):
    """
    Analyses Resource
    """
    program=None
    programversion=None
    name = None
    sourcename = None
    description = None
    algorithm = None
    sourceversion = None
    sourceuri = None
    date_executed = None

    def get(self, id=None):
        print(f'GET {request.path}\nGetting analyses {id}')
        self._check_data(request.args)
        from biobarcoding.services.analyses import read_analyses
        response, code = read_analyses(
            analysis_id=id,
            name=self.name,
            program=self.program,
            programversion=self.programversion,
            algorithm=self.algorithm,
            sourcename=self.sourcename,
            sourceversion=self.sourceversion,
            sourceuri=self.sourceuri,
            description=self.description)
        return make_response(jsonify(response), code)


    def post(self):
        print(f'POST {request.path}\nCreating analyses')
        self._check_data(request.json)
        from biobarcoding.services.analyses import create_analyses
        response, code = create_analyses(
            program=self.program,
            programversion=self.programversion,
            name=self.name,
            sourcename=self.sourcename,
            description=self.description,
            algorithm=self.algorithm,
            sourceversion=self.sourceversion,
            sourceuri=self.sourceuri,
            date_executed=self.date_executed)
        return make_response(jsonify(response), code)


    def put(self, id):
        print(f'POST {request.path}\nCreating analyses')
        self._check_data(request.json)
        from biobarcoding.services.analyses import update_analyses
        response, code = update_analyses(id,
            program=self.program,
            programversion=self.programversion,
            name=self.name,
            sourcename=self.sourcename,
            description=self.description,
            algorithm=self.algorithm,
            sourceversion=self.sourceversion,
            sourceuri=self.sourceuri,
            date_executed=self.date_executed)
        return make_response(jsonify(response), code)


    def delete(self, id):
        print(f'DELETE {request.path}\nDeleting analyses {id}')
        self._check_data(request.args)
        from biobarcoding.services.analyses import delete_analyses
        response, code = delete_analyses(id,
            name=self.name,
            program=self.program,
            programversion=self.programversion,
            algorithm=self.algorithm,
            sourcename=self.sourcename,
            sourceversion=self.sourceversion,
            sourceuri=self.sourceuri,
            description=self.description)
        return make_response(jsonify(response), code)


    def _check_data(self, data):
        if data:
            if 'program' in data and data['program']:
                self.program = data['program']
            if 'programversion' in data and data['programversion']:
                self.programversion = data['programversion']
            if 'name' in data and data['name']:
                self.name = data['name']
            if 'sourcename' in data and data['sourcename']:
                self.sourcename = data['sourcename']
            if 'description' in data and data['description']:
                self.description = data['description']
            if 'algorithm' in data and data['algorithm']:
                self.algorithm = data['algorithm']
            if 'sourceversion' in data and data['sourceversion']:
                self.sourceversion = data['sourceversion']
            if 'sourceuri' in data and data['sourceuri']:
                self.sourceuri = data['sourceuri']
            if 'date_executed' in data and data['date_executed']:
                self.date_executed = data['date_executed']
        print(f'DATA: {data}')


analyses_view = AnalysesAPI.as_view('api_analyses')
bp_analyses.add_url_rule(
    bcs_api_base + '/analyses/',
    view_func=analyses_view,
    methods=['GET','POST']
)
bp_analyses.add_url_rule(
    bcs_api_base + '/analyses/<int:id>',
    view_func=analyses_view,
    methods=['GET','PUT','DELETE']
)
