import os.path
import time

from Bio import Phylo
# from Bio.Phylo.BaseTree import Tree, Clade
from dendropy import Tree, TreeList, Node, Taxon

from .sequences import Service as SeqService
from .analyses import Service as AnsisService
from ..meta.ontologies import get_type_id
from ... import get_or_create, log_exception, get_orm_params
from ....common.helpers import zip_files
from ....db_models.bioinformatics import PhylogeneticTree
from ....db_models.chado import Phylotree, Phylonode, Feature, Phylonodeprop


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
                _ = [node.label or node.taxon.label if node.taxon else node.label for node in tree]
                feats = dict(self.db.query(Feature.uniquename, Feature.feature_id)
                             .distinct(Feature.feature_id).filter(Feature.uniquename.in_(_)).all())
                self.db.add_all(self.tree2phylonodes(pt.phylotree_id, tree.seed_node, None, [0], feats))
                # phylonodes = self.tree2phylonodes(pt.phylotree_id, tree.root, None, [0])
        except Exception as e:
            log_exception(e)
            raise Exception(f'IMPORT phylotress: The file {os.path.basename(infile)} could not be imported. (unmanageable)')
        return content, count

    def tree2phylonodes(self, phylotree_id, node, parent_id=None, index=[0], features={}):
        new_entries = []
        n_name = node.label or node.taxon.label if node.taxon else node.label
        # n_name = node.name or node.id
        phylonode = Phylonode(phylotree_id=phylotree_id,
                              parent_phylonode_id=parent_id,
                              feature_id=features.get(n_name, None),
                              label=n_name,
                              # label=node.id or node.name,
                              distance=node.edge_length,
                              # distance=node.branch_length,
                              left_idx=index[0],
                              right_idx=index[0] + 1)
        for k in ('height', 'posterior', 'height_95%_HPD'):
            v = node.annotations.get_value(k)
            if v:
                if not phylonode.phylonode_id:
                    self.db.add(phylonode)
                    self.db.flush()
                new_entries.append(Phylonodeprop(phylonode_id=phylonode.phylonode_id,
                                                 type_id=1, value=f'{k}:{v}'))    # TODO: tag the values
        # Check for children
        index[0] += 1
        _ = node.child_nodes()
        # _ = node.clades
        if len(_) > 0:
            if not phylonode.phylonode_id:
                self.db.add(phylonode)
                self.db.flush()
            for clade in _:
                new_entries += self.tree2phylonodes(phylotree_id, clade, phylonode.phylonode_id, index, features)
            phylonode.right_idx = index[0]
        index[0] += 1
        return [phylonode] + new_entries

    ##
    # EXPORT
    ##

    def data2file(self, data: list, outfile, format: str, values={}, **kwargs) -> int:
        files, count = [], 0
        for ans in data:
            trees = TreeList()
            _ = kwargs.get('phylotree_id')
            pts = ans.phylotrees if not _ \
                else [t for t in ans.phylotrees if t.phylotree_id in _ or t.phylotree_id == _]
            for pt in pts:
                root = self.db.query(Phylonode).filter(Phylonode.phylotree_id == pt.phylotree_id,
                                                       Phylonode.parent_phylonode_id.is_(None)).one()
                # rooted = bool(root.feature_id or root.label)
                # trees.append(Tree(root=self.chado2biopy(root, values.get('header')), rooted=rooted, id=pt.name, name=pt.name))
                trees.append(Tree(seed_node=self.chado2biopy(root, values.get('header')), label=pt.name))
                count += 1
            _file = outfile if len(data) < 2 else f'{outfile}_{ans.name}_{ans.analysis_id}'
            try:
                trees.write(path=_file, schema=format)
                # Phylo.write(trees, _file, format)
            except Exception as e:
                if format not in self.formats:
                    raise Exception('The specified format is not available.\n'
                                    + 'The supported formats for phylogenetic trees are: ' + str(self.formats))
                if len(trees) > 1:
                    _fs = []
                    for t in trees:
                        _f = f'{_file}_{t.label}'
                        # _f = f'{_file}_{t.name}_{t.phylotree_id}'
                        Phylo.write(t, _f, format)
                        _fs.append(_f)
                    zip_files(_file, _fs)
                else:
                    raise e
            files.append(_file)
        if len(files) > 1:
            zip_files(outfile, files)
        return count

    def chado2biopy(self, node, header: str = None) -> Node:
        _label = node.label
        if header and node.feature_id:
            # retrieve and load the label
            _ = self.db.query(Feature).filter(Feature.feature_id == node.feature_id).one()
            _label = SeqService().seqs_header_parser([_], header).get(_label, _label)
        clade = Node(edge_length=node.distance, label=_label or '', taxon=Taxon(label=_label or ''))
        _ann = [_.value.split(':', 1) for _ in self.db.query(Phylonodeprop)
            .filter(Phylonodeprop.phylonode_id == node.phylonode_id).all()]
        if _ann:
            [clade.annotations.add_new(_[0], _[-1]) for _ in _ann]
        # clade = Clade(branch_length=node.distance, name=_label or '')
        children = self.db.query(Phylonode).filter(Phylonode.parent_phylonode_id == node.phylonode_id)
        clade.set_child_nodes([self.chado2biopy(n, header) for n in children.all()])
        # clade.clades = [self.chado2biopy(n, header) for n in children.all()]
        return clade

    ##
    # GETTER AND OTHERS
    ##

    def aux_filter(self, _filter: dict) -> list:
        clauses = []
        from ....rest import filter_parse

        if _filter.get('feature_id'):
            _ids = self.db.query(Phylotree.analysis_id).join(Phylonode) \
                .filter(filter_parse(Phylonode, [{'feature_id': _filter.get('feature_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if _filter.get('organism_id'):
            from ....db_models.chado import Feature
            _ids = self.db.query(Phylotree.analysis_id).join(Phylonode).join(Feature) \
                .filter(filter_parse(Feature, [{'organism_id': _filter.get('organism_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if _filter.get("prop_cvterm_id"):
            from ....db_models.chado import Phylotreeprop
            _ids = self.db.query(Phylotree.analysis_id).join(Phylotreeprop) \
                .filter(filter_parse(Phylotreeprop, [{'type_id': _filter.get('prop_cvterm_id')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        from ....db_models.chado import Analysis
        if _filter.get("program"):
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'program': _filter.get('program')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if _filter.get("programversion"):
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'programversion': _filter.get('programversion')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        if _filter.get("algorithm"):
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, [{'algorithm': _filter.get('algorithm')}]))
            clauses.append(self.orm.analysis_id.in_(_ids))

        from datetime import datetime
        if _filter.get("added-from"):
            _filter["added-from"]['unary'] = datetime.strptime(_filter.get("added-from")['unary'], '%Y-%m-%d')
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, {'timeexecuted': _filter.get("added-from")}))
            clauses.append(self.orm.analysis_id.in_(_ids))
        if _filter.get("added-to"):
            _filter["added-to"]['unary'] = datetime.strptime(_filter.get("added-to")['unary'], '%Y-%m-%d')
            _ids = self.db.query(Analysis.analysis_id) \
                .filter(filter_parse(Analysis, {'timeexecuted': _filter.get("added-to")}))
            clauses.append(self.orm.analysis_id.in_(_ids))

        return clauses + super(Service, self).aux_filter(_filter)
