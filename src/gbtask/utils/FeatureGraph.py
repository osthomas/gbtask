from typing import Callable, Hashable, Self, Set

import networkx as nx
from Bio.SeqFeature import SeqFeature


class FeatureGraph:
    """
    Mixin class for utility functions to handle SeqFeature objects associated
    with a nx.Graph node.
    """

    def __init__(self, graph: nx.Graph | nx.DiGraph):
        self.graph: nx.Graph | nx.DiGraph = graph

    @property
    def features(self):
        """
        Return a list of features of the underlying feature graph.
        """
        return [self.graph.nodes[n]["feature"] for n in self.graph]

    @property
    def node_groups(self) -> list[set[Hashable]]:
        """
        Return a list of sets of node keys from the feature graph.
        Each set contains the node keys of one connected component in the graph,
        ie., nodes for groups of features.
        To access the nodes themselves, use `self.graph[node_id]`.

        See also
        --------
        :meth:`feature_groups`
        """
        for fun in (nx.connected_components, nx.weakly_connected_components):
            try:
                return list(fun(self.graph))
            except nx.NetworkXNotImplemented:
                pass
        raise NotImplementedError("Expected nx.Graph or nx.DiGraph as self.graph")

    @property
    def feature_groups(self) -> list[list[SeqFeature]]:
        """
        Return a list of lists of features from the feature graph.
        Each sublist contains the features of one connected component in the
        graph, ie., groups of features.

        See also
        --------
        :meth:`node_groups`
        """
        feature_groups = []
        for nodes in self.node_groups:
            feature_group = []
            for node in nodes:
                feature_group.append(self.graph.nodes[node]["feature"])
            feature_groups.append(feature_group)
        return feature_groups

    def add_feature_node(self, feature: SeqFeature):
        """
        Add a new feature node to feature graph.
        A SeqFeature is not hashable and cannot be added to the graph directly.
        """
        node_id = hex(id(feature))
        new_id = node_id
        while new_id in self.graph:
            # Make absolutely sure the ID is unique
            # NOTE: this means the same feature can be added to the graph
            # multiple times!
            increment = 0
            new_id = node_id + str(increment)
        self.graph.add_node(new_id, feature=feature)
        return new_id

    def per_group(self, callback: Callable[[Self, Set[Hashable]], None]):
        """
        Execute `callback` on all nodes groups.

        `callback` receives 2 arguments:
            1. a reference to the "FeatureGraph" instance
            2. the set of node IDs for one feature group
        """
        for nodes in self.node_groups:
            callback(self, nodes)
