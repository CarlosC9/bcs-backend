import os.path
import time

from . import BosService
from ..meta.ontologies import get_type_id
from ...main import get_orm
from ....db_models import DBSessionChado, DBSession
from ....db_models.bioinformatics import PhylogeneticTree
from ....db_models.chado import Phylonode, Feature
from ... import get_bioformat, get_or_create, log_exception


##
# PHYLOTREE SERVICE
##
class Service(BosService):

    def __init__(self):
        super(Service, self).__init__()
        self.db = DBSessionChado
        self.orm = get_orm('phylotrees')
        self.bos = 'phylogenetic-tree'

    ##
    # CREATE
    ##

    def prepare_values(self, **values):

        if not values.get('type_id') and values.get('type'):
            try:
                values['type_id'] = get_type_id(type=values.get('type'), subtype=values.get('subtype'))
            except:
                pass

        # Analysis row could exist for jobs, so get or create
        ansis_trigger = ('job_id', 'program', 'programversion', 'sourcename')
        if not values.get('analysis_id') and any(key in values.keys() for key in ansis_trigger):
            from ..meta.analyses import Service as AnsisService
            try:
                unique_keys = ['job_id'] if values.get('job_id') else ('program', 'programversion', 'sourcename')
                content, count = AnsisService.get_query(**{k:values[k] for k in unique_keys if k in values})
                content = content.one()
            except:
                t_name = values.pop('name')   # Avoid to spread the treename to analysis
                content, count = AnsisService.create(**values)
                values['name'] = t_name
            values['analysis_id'] = content.analysis_id

        return super(Service, self).prepare_values(**values)

    def check_values(self, **values):

        if not values.get('dbxref_id'):
            from ....db_models.chado import Db, Dbxref
            values['dbxref_id'] = get_or_create(self.db, Dbxref,
                                                db_id=self.db.query(Db).filter(Db.name == 'null').one().db_id,
                                                accession=f'phylotree:{values.get("name")}',
                                                version=time.strftime("%Y %b %d %H:%M:%S")).dbxref_id

        if not values.get('type_id'):
            try:
                values['type_id'] = get_type_id(type='phylotree')
            except:
                pass

        return super(Service, self).check_values(**values)

    def after_create(self, new_object, **values):
        # tree to bcs
        phylo = get_or_create(DBSession, PhylogeneticTree,
                              native_id=new_object.phylotree_id,
                              native_table='phylotree',
                              name=new_object.name)
        return values

    ##
    # DELETE
    ##

    def after_delete(self, *content, **kwargs):
        ids = [t.phylotree_id for t in content]
        query = DBSession.query(PhylogeneticTree).filter(PhylogeneticTree.native_id.in_(ids))
        return query.count()
        # TODO: check why all rows are deleted in bcs without filtering
        # return query.delete(synchronize_session='fetch')

    ##
    # IMPORT
    ##

    def import_file(self, input_file, format=None, **kwargs):
        ## TODO: NGD newick phylotree import
        # phylotree: store the input_file content into a phylotreeprop EDAM_newick ?
        # phylonode: type_id ? (root,leaf,internal)
        content, count = None, 0
        format = get_bioformat(input_file, format)
        try:
            # Check phylotree file format
            from Bio import Phylo
            if format == 'nexus':
                from dendropy import TreeList
                from io import StringIO
                trees = TreeList.get_from_path(input_file, schema="nexus")
                tree = Phylo.read(StringIO(trees[-1].as_string("nexus")), "nexus")
            else:
                # TODO: try formats ('newick', 'nexus', 'nexml', 'phyloxml', 'cdao') ?
                tree = Phylo.read(input_file, format)
        except:
            raise Exception(f'IMPORT phylotress: The file {os.path.basename(input_file)} could not be imported. (unknown format)')
        try:
            if not kwargs.get('name'):
                kwargs['name'] = os.path.basename(input_file)
            # Create the new phylotree
            content, count = self.create(**kwargs)
            # Get phylonodes insertion
            phylonodes = self.tree2phylonodes(content.phylotree_id, tree.root, None, [0])
            count += 1
        except Exception as e:
            log_exception(e)
            raise Exception(f'IMPORT phylotress: The file {os.path.basename(input_file)} could not be imported. (unmanageable)')
        return content, count


    def tree2phylonodes(self, phylotree_id, node, parent_id=None, index=[0]):
        phylonodes = []
        f = self.db.query(Feature).filter(Feature.uniquename == node.name).first()
        try:
            feature_id = f.feature_id
        except:
            feature_id = None
        phylonode = get_or_create(self.db, Phylonode,
                                  phylotree_id=phylotree_id,
                                  parent_phylonode_id=parent_id,
                                  feature_id=feature_id,
                                  label=node.name,
                                  distance=node.branch_length,
                                  left_idx=index[0],
                                  right_idx=index[0] + 1)
        # Check for children
        index[0] += 1
        if len(node.clades) > 0:
            for clade in node.clades:
                phylonodes += self.tree2phylonodes(phylotree_id, clade, phylonode.phylonode_id, index)
            phylonode.right_idx = index[0]
            self.db.merge(phylonode)
        index[0] += 1
        return [phylonode] + phylonodes


    ##
    # EXPORT
    ##

    def export_file(self, format='newick', **kwargs):
        content, count = f'/tmp/output_ngd.{format}', 0
        # NGD newick phylotree export
        try:
            t_id = self.get_query(purpose='share', **kwargs)[0].one().phylotree_id
            self.tree2file(t_id, format, content)
            count += 1
        except Exception as e:
            log_exception(e)
            raise Exception('EXPORT phylotrees: The phylotree could not be exported.')
        return content, count


    def tree2file(self, phylotree_id, format, output_file):
        format = get_bioformat(output_file, format)
        # TODO: build Bio.Phylo.BaseTree from chado, and Phylo.write(tree, output_file, format)
        if format:
            # dump to newick
            root = self.db.query(Phylonode).filter(Phylonode.phylotree_id == phylotree_id,
                                                         Phylonode.parent_phylonode_id == None).one()
            result = self.tree2newick(root)
            with open("%s.newick" % output_file, "w") as file:
                file.write(result)
            # convert from newick to another format
            if format in ["nexus", "nexml"] :
                from dendropy import Tree
                tree = Tree.get(path="%s.newick" % output_file, schema="newick")
                tree.write(path=output_file, schema=format)
            elif format != "newick":
                from Bio import Phylo
                Phylo.convert("%s.newick" % output_file, "newick", output_file, format)
            else:
                os.rename("%s.newick" % output_file, output_file)
        return output_file

    def tree2newick(self, node, is_root=True):
        result = ''
        children = self.db.query(Phylonode).filter(Phylonode.parent_phylonode_id == node.phylonode_id)
        if children.count() > 0:
            result += '(' if is_root else '('
            i = 1
            for n in children.all():
                result += self.tree2newick(n, False)
                result += ',' if i < children.count() else ')'
                i += 1
        if not is_root and node.label or node.distance:
            label = node.label if node.label else ''
            distance = node.distance if node.distance else '0.00000'
            result += f'{label}:{distance}'
        result += ';' if is_root else ''
        return result

    ##
    # GETTER AND OTHERS
    ##

    def aux_filter(self, filter):
        clauses = []
        from ....rest import filter_parse

        if 'feature_id' in filter:
            _ids = self.db.query(Phylonode.phylotree_id) \
                .filter(filter_parse(Phylonode, [{'feature_id': filter.get('feature_id')}]))
            clauses.append(self.orm.phylotree_id.in_(_ids))

        if 'organism_id' in filter:
            from ....db_models.chado import Feature
            _ids = self.db.query(Feature.feature_id) \
                .filter(filter_parse(Feature, [{'organism_id': filter.get('organism_id')}]))
            _ids = self.db.query(Phylonode.phylotree_id) \
                .filter(Phylonode.feature_id.in_(_ids))
            clauses.append(self.orm.phylotree_id.in_(_ids))

        if "prop_cvterm_id" in filter:
            from ....db_models.chado import Phylotreeprop
            _ids = self.db.query(Phylotreeprop.phylotree_id) \
                .filter(filter_parse(Phylotreeprop, [{'type_id': filter.get('prop_cvterm_id')}]))
            clauses.append(self.orm.phylotree_id.in_(_ids))

        from ....db_models.chado import Analysis
        if "program" in filter:
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'program': filter.get('program')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if "programversion" in filter:
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'programversion': filter.get('programversion')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if "algorithm" in filter:
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'algorithm': filter.get('algorithm')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        from datetime import datetime
        if "added-from" in filter:
            filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, {'timeexecuted': filter.get("added-from")}))
            clauses.append(self.orm.analysis_id.in_(_ids))
        if "added-to" in filter:
            filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, {'timeexecuted': filter.get("added-to")}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        return clauses
