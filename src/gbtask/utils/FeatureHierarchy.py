import itertools
import logging
from typing import Collection, Hashable, Set

import networkx as nx
from Bio.SeqFeature import SeqFeature
from gbtask.utils import futils
from gbtask.utils.FeatureGraph import FeatureGraph
from gbtask.utils.FeatureIDManager import FeatureIDManager
from gbtask.utils.TypeTree import DefaultTypeTree, TypeTree

log = logging.getLogger(__name__)


class FeatureHierarchy(FeatureGraph):
    def __init__(
        self,
        graph: nx.Graph | nx.DiGraph,
        add_root: bool = False,
        type_tree: TypeTree = DefaultTypeTree,
    ):
        """
        Initialize a FeatureHierarchy from a pre-connected feature graph.


        Parameters
        ----------

        feature_graph
            Pre-connected feature graph.
            Nodes in `feature_graph` must already be meaningfully connected.
            FeatureHierarchy does not perform grouping, just re-organization of
            edges within groups to reflect the underlying type hierarchy.

        add_root
            If True, root type features are automatically added to groups in
            `feature_graph` if necessary.
            See :meth:`_add_root_node_to` for details.

        type_tree
            A Tree representing the hierarchy of feature types.
        """
        # possibly undirected before hierarchization, directed after
        # Conversion nx.Graph -> nx.Digraph replaces each undirected edge by
        # two directional edges
        super().__init__(nx.DiGraph(graph))
        self._type_tree: TypeTree = DefaultTypeTree
        self._hierarchize(add_root)

    def _is_parent(self, parent: SeqFeature, child: SeqFeature):
        """
        Return True if `parent` is a parent feature type of `child` according
        to the underlying feature type tree.
        """
        if parent.type not in self._type_tree._tree.nodes:
            log.info(f"{parent.type} is not in the feature hierarchy.")
        return parent.type in self._type_tree.parent_types(child.type)

    def _get_edge(self, node1: Hashable, node2: Hashable):
        """
        Return the directed edge between `node1` and `node2` based on the
        underlying feature type hierarchy.
        """
        parent_feature = self.graph.nodes[node1]["feature"]
        child_feature = self.graph.nodes[node2]["feature"]
        if self._is_parent(parent_feature, child_feature):
            return (node1, node2)
        elif self._is_parent(child_feature, parent_feature):
            return (node2, node1)
        else:
            return tuple()

    def _hierarchize_group(self, nodes: Collection[Hashable]):
        """
        Hierarchize one connected component of the complete feature graph.
        """
        # start from a clean slate: remove all edges involving the nodes
        edges = list(self.graph.edges(nbunch=nodes))
        self.graph.remove_edges_from(edges)
        # Separate nodes into root and nonroot nodes
        root_nodes = set()
        nonroot_nodes = set()
        for node in nodes:
            if self.graph.nodes[node]["feature"].type == self._type_tree.root_type:
                root_nodes.add(node)
            else:
                nonroot_nodes.add(node)
        log.debug(f"nodes: {nodes} -- roots: {root_nodes} -- nonroots: {nonroot_nodes}")
        for parent_node, child_node in itertools.combinations(nodes, 2):
            _debug = f"Checking link {parent_node} <-> {child_node}"
            if parent_node == child_node:
                _debug += " - skipped"
                continue
            edge = self._get_edge(parent_node, child_node)
            if edge:
                _debug += f" - edge {edge}"
                self.graph.add_edge(*edge)
            log.debug(_debug)
        # A node without a parent at this point is either:
        # * a root feature, in which case we do nothing
        # * or it is a dangling feature that is not part of the defined type
        #   hierarchy, so we link it directly to the root feature,
        # * or it is a feature that could have a parent, but none is
        #   available in this tree, so we also link it directly to the root
        #   feature (eg. an 'exon' feature in a group without an 'mRNA')
        noparent_nodes = [n for n in nonroot_nodes if self.graph.in_degree(n) == 0]
        log.debug(f"starting rescue for: {noparent_nodes}")
        for root in root_nodes:
            for noparent in noparent_nodes:
                _debug = f"Checking link {root} -> {noparent}"
                if root != noparent:
                    edge = (root, noparent)
                    self.graph.add_edge(root, noparent)
                    _debug += f" - add edge {edge}"
                    log.debug(_debug)
        return nodes

    def _hierarchize(self, add_root: bool = False):
        """
        Hierarchize the feature graph

        Modifies `self.graph` in place.
        """
        for nodes in self.node_groups:
            if add_root and self.n_of_type(self._type_tree.root_type, nodes) == 0:
                new_node = self._add_root_node_to(nodes)
                if new_node is not None:
                    nodes.add(new_node)
            self._hierarchize_group(nodes)

    def n_of_type(self, type: str, nodes: Set[Hashable]):
        """
        Return the number of features of `type` in the set of `nodes`.
        """
        n = 0
        for node in nodes:
            n += self.graph.nodes[node]["feature"] == type
        return n

    def _add_root_node_to(self, nodes: Set[Hashable]):
        """
        Add a root feature to a feature group without one.
        The created feature will cover the complete range of all other features
        in the group that are part of the type tree, including gaps.
        A root feature is only added if the group contains at least one feature
        that is associated with the root feature in the underlying type tree.
        This avoids creation of spurious root features for isolated, dangling
        features in the feature graph (eg. stray primer_bind features etc.)

        Returns
        -------
        Key of the newly created node, or None if no root node was created.
        """
        features = []
        for n in nodes:
            feature = self.graph.nodes[n]["feature"]
            if feature.type == self._type_tree.root_type:
                # This is a root type, don't need to add one
                return None
            if feature.type in self._type_tree._tree.nodes:
                features.append(feature)
        if features:
            new_root = futils.create_spanning(self._type_tree.root_type, features)
            nid = len(self.graph)
            # TODO: come up with a better way to handle node IDs
            while nid in self.graph:
                nid += 1
            self.graph.add_node(nid, feature=new_root)
            return nid
        else:
            return None

    def set_feature_ids(self):
        """
        Set feature IDs for features in the hierarchy.

        Features IDs are set based on the features' types with a running count.
        Two classes are distinguished:
            1. Top-level features at the root level of the feature hierarchy (ie.,
               these features have no parent features).
               They receive an ID based on their feature type, numbered
               sequentially according to start position.
            2. Non-top-level features which have parent features.
               They are numbered within the group and then prefixed with the ID
               of their group's top-level feature(s).
        """
        # Set IDs for top-level features without parents first
        top_nodes = [n for n, deg in self.graph.in_degree() if deg == 0]
        top_features = [self.graph.nodes[n]["feature"] for n in top_nodes]
        top_mgr = FeatureIDManager(top_features)
        top_mgr.set_feature_ids()
        for group in self.node_groups:
            top_group = []  # top level nodes in the group (no parents)
            other_group = []  # other features in the group (have parents)
            for n in group:
                if self.graph.in_degree(n) == 0:
                    top_group.append(n)
                else:
                    other_group.append(n)
            # Construct group prefix from all root features in this group
            # If there is more than one root feature, the prefix is built from
            # concatenating individual root feature IDs
            group_prefix = [self.graph.nodes[n]["feature"].id for n in top_group]
            group_prefix = "_".join(sorted(group_prefix))  # sort for determinism
            group_features = [self.graph.nodes[n]["feature"] for n in other_group]
            group_mgr = FeatureIDManager(group_features, prefix=group_prefix)
            group_mgr.set_feature_ids()

    def _set_gff_qualifiers(self):
        """
        Set Parents for the GFF format.
        GFF is based on ID= and Parent= qualifiers.
        A feature can have multiple parents.
        This information is readily available from the hierarchy graph, but not
        from a flat feature list, eg. after export to Genbank.
        Bake the parent information into the feature qualifiers.
        """
        for node in self.graph:
            feature = self.graph.nodes[node]["feature"]
            if "Parent" not in feature.qualifiers:
                feature.qualifiers["Parent"] = []
            # sorted because ordering of predecessors is non-deterministic
            for parent_node in sorted(self.graph.predecessors(node)):
                parent_feature = self.graph.nodes[parent_node]["feature"].id
                if parent_feature not in feature.qualifiers["Parent"]:
                    feature.qualifiers["Parent"].append(parent_feature)

    def _set_gtf_qualifiers(self):
        """
        Set grouping qualifiers for the GTF format.
        GTF is based on gene_id and transcript_id qualifiers.
        The hierarchy graph must contain annotations to identify which feature
        types should be taken to define gene_id and transcript_id qualifiers.
        gene_id should be set to the ID of the parent gene feature.
        transcript_id should be set to the ID of a parent RNA feature.
        """
        for node, data in self.graph.nodes(data=True):
            feature = data["feature"]
            for gtf_level in self._type_tree.gtf_levels:
                types = self._type_tree.gtf_level_types(gtf_level)
                if feature.type in types:
                    for successor in nx.dfs_tree(self.graph, node).nodes:
                        child = self.graph.nodes[successor]["feature"]
                        key = f"{gtf_level}"
                        if key not in child.qualifiers:
                            child.qualifiers[key] = []
                        if feature.id not in child.qualifiers[key]:
                            child.qualifiers[key].append(feature.id)

    def to_features(self):
        """
        Augment feature qualifiers with hierarchical information.
        """
        self.set_feature_ids()
        self._set_gff_qualifiers()
        self._set_gtf_qualifiers()
