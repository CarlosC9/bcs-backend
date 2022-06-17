import os.path
import time

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade

from .analyses import Service as AnsisService
from ..meta.ontologies import get_type_id
from ... import get_or_create, log_exception, get_orm_params
from ....common.helpers import zip_files
from ....db_models.bioinformatics import PhylogeneticTree
from ....db_models.chado import Phylotree, Phylonode, Feature


##
# PHYLOTREE SERVICE
##
class Service(AnsisService):

    def __init__(self):
        super(Service, self).__init__()
        self.obj_type = 'phylogenetic-tree'
        self.fos = PhylogeneticTree
        self.formats = ['fasta', 'nexus', 'nexml', 'phylip', 'phyloxml', 'newick']     # newick at the end because it swallows everything

    ##
    # CREATE
    ##

    def prepare_external_values(self, **values) -> dict:

        if not values.get('type_id') and values.get('type'):
            try:
                values['type_id'] = get_type_id(type=values.get('type'), subtype=values.get('subtype'))
            except:
                pass

        if not values.get('type_id'):
            try:
                values['type_id'] = get_type_id(type='phylotree')
            except:
                pass

        if not values.get('dbxref_id'):
            from ....db_models.chado import Db, Dbxref
            values['dbxref_id'] = get_or_create(self.db, Dbxref,
                                                db_id=self.db.query(Db).filter(Db.name == 'null').one().db_id,
                                                accession=f'phylotree:{values.get("name")}',
                                                version=time.strftime("%Y %b %d %H:%M:%S")).dbxref_id

        return super(Service, self).prepare_external_values(**values)

    ##
    # IMPORT
    ##

    def read_infile(self, file, _format) -> any:
        from dendropy import Tree
        _ty = Tree.yield_from_files(files=[file], schema=_format)
        for _ in _ty:
            return Tree.yield_from_files(files=[file], schema=_format)
            # return Phylo.parse(file, _format)

    def import_file(self, infile, format=None, **kwargs):
        # TODO ? phylonode:type_id (root,leaf,internal) ?
        content_file, _format = self.check_infile(infile, format)
        try:
            # Set missing default values
            kwargs['programversion'] = kwargs.get('programversion') or '(Imported file)'
            kwargs['sourcename'] = kwargs.get('sourcename') or os.path.basename(infile)
            if not kwargs.get('name'):
                kwargs['name'] = os.path.basename(infile)
            # Create the new phylotree
            content, count = self.create(**kwargs)
            _name = kwargs.pop('name')
            for tree in content_file:
                pt = get_or_create(self.db, Phylotree, analysis_id=content.analysis_id,
                                   name=tree.label or _name,
                                   # name=tree.id or tree.name or _name,
                                   **get_orm_params(Phylotree, **self.prepare_external_values(**kwargs)))
                # Get phylonodes insertion
                phylonodes = self.tree2phylonodes(pt.phylotree_id, tree.seed_node, None, [0])
                # phylonodes = self.tree2phylonodes(pt.phylotree_id, tree.root, None, [0])
        except Exception as e:
            log_exception(e)
            raise Exception(f'IMPORT phylotress: The file {os.path.basename(infile)} could not be imported. (unmanageable)')
        return content, count

    def tree2phylonodes(self, phylotree_id, node, parent_id=None, index=[0]):
        phylonodes = []
        try:
            n_name = node.label or node.taxon.label
            # n_name = node.name or node.id
            feature_id = self.db.query(Feature).filter(
                Feature.uniquename == n_name).one().feature_id
        except:
            n_name = feature_id = None
        phylonode = get_or_create(self.db, Phylonode,
                                  phylotree_id=phylotree_id,
                                  parent_phylonode_id=parent_id,
                                  feature_id=feature_id,
                                  label=n_name,
                                  # label=node.id or node.name,
                                  distance=node.edge_length,
                                  # distance=node.branch_length,
                                  left_idx=index[0],
                                  right_idx=index[0] + 1)
        # Check for children
        index[0] += 1
        _ = node.child_nodes()
        # _ = node.clades
        if len(_) > 0:
            for clade in _:
                phylonodes += self.tree2phylonodes(phylotree_id, clade, phylonode.phylonode_id, index)
            phylonode.right_idx = index[0]
            self.db.merge(phylonode)
        index[0] += 1
        return [phylonode] + phylonodes

    ##
    # EXPORT
    ##

    def data2file(self, _phys: list, outfile, format: str, **kwargs) -> int:
        trees, count = [], 0
        files = []
        for ans in _phys:
            _ = kwargs.get('phylotree_id')
            pts = ans.phylotrees if not _ \
                else [t for t in ans.phylotrees if t.phylotree_id in _ or t.phylotree_id == _]
            for pt in pts:
                root = self.db.query(Phylonode).filter(Phylonode.phylotree_id == pt.phylotree_id,
                                                       Phylonode.parent_phylonode_id.is_(None)).one()
                rooted = bool(root.feature_id or root.label)
                trees.append(Tree(root=self.tree2biopy(root), rooted=rooted, id=pt.name, name=pt.name))
                count += 1
            _file = outfile if len(_phys) < 2 else f'{ans.name or outfile}_{ans.analysis_id}'
            try:
                Phylo.write(trees, _file, format)
                # if format not in ["nexus", "nexml"]:
                #     Phylo.write(trees, "%s.newick" % outfile, 'newick')
                #     from dendropy import Tree as DenTree
                #     tree = DenTree.get(path="%s.newick" % outfile, schema="newick")
                #     tree.write(path=file, schema=format)
            except Exception as e:
                if format not in self.formats:
                    raise Exception('The specified format is not available.\n'
                                    + 'The supported formats for phylogenetic trees are: ' + str(self.formats))
                if len(trees) > 1:
                    _fs = []
                    for t in trees:
                        _f = f'{t.name or _file}_{t.phylotree_id}'
                        Phylo.write(t, _f, format)
                        _fs.append(_f)
                    zip_files(_file, _fs)
                else:
                    raise e
            files.append(_file)
        if len(files) > 1:
            zip_files(outfile, files)
        return count

    def tree2biopy(self, node, label_type=None) -> Clade:
        if label_type:
            # TODO: retrieve and change the label
            pass
        clade = Clade(branch_length=node.distance, name=node.label or '')
        children = self.db.query(Phylonode).filter(Phylonode.parent_phylonode_id == node.phylonode_id)
        clade.clades = [self.tree2biopy(n) for n in children.all()]
        return clade

    ##
    # GETTER AND OTHERS
    ##

    def aux_filter(self, filter):
        clauses = []
        from ....rest import filter_parse

        if filter.get('feature_id'):
            _ids = self.db.query(Phylotree.analysis_id).join(Phylonode) \
                .filter(filter_parse(Phylonode, [{'feature_id': filter.get('feature_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get('organism_id'):
            from ....db_models.chado import Feature
            _ids = self.db.query(Phylotree.analysis_id).join(Phylonode).join(Feature) \
                .filter(filter_parse(Feature, [{'organism_id': filter.get('organism_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get("prop_cvterm_id"):
            from ....db_models.chado import Phylotreeprop
            _ids = self.db.query(Phylotree.analysis_id).join(Phylotreeprop) \
                .filter(filter_parse(Phylotreeprop, [{'type_id': filter.get('prop_cvterm_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        from ....db_models.chado import Analysis
        if filter.get("program"):
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'program': filter.get('program')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get("programversion"):
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'programversion': filter.get('programversion')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if filter.get("algorithm"):
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'algorithm': filter.get('algorithm')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        from datetime import datetime
        if filter.get("added-from"):
            filter["added-from"]['unary'] = datetime.strptime(filter.get("added-from")['unary'], '%Y-%m-%d')
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, {'timeexecuted': filter.get("added-from")}))
            clauses.append(self.orm.analysis_id.in_(_ids))
        if filter.get("added-to"):
            filter["added-to"]['unary'] = datetime.strptime(filter.get("added-to")['unary'], '%Y-%m-%d')
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, {'timeexecuted': filter.get("added-to")}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(filter)
