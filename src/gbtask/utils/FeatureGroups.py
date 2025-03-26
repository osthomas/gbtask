import itertools
import logging
from typing import Sequence

import networkx as nx
from Bio.SeqFeature import SeqFeature
from gbtask.utils import futils
from gbtask.utils.FeatureGraph import FeatureGraph

log = logging.getLogger(__name__)


class FeatureGroups(FeatureGraph):
    """
    Arrange a collection of `SeqFeature`s into a graph structure based on their
    connectivity.
    """

    def __init__(
        self,
        features: list[SeqFeature],
        node_label_from_id: bool = False,
    ):
        """
        Initialize FeatureGroups from a list of features.

        Parameters
        ----------

        features
            list of SeqFeatures in the feature graph

        node_label_from_id
            If True, assign node labels based on the feature id to facilitate
            access. Feature IDs must be unique. If False, nodes are labelled
            sequentially.

        """
        self.graph: nx.Graph = nx.Graph()
        for i, feature in enumerate(features):
            if node_label_from_id:
                if feature.id not in self.graph.nodes:
                    self.graph.add_node(feature.id, feature=feature)
                else:
                    raise ValueError(f"Feature ID {feature.id} is not unique")
            else:
                self.add_feature_node(feature)
        super().__init__(self.graph)

    @classmethod
    def _connectable_by_qualifier(
        cls,
        feature1: SeqFeature,
        feature2: SeqFeature,
        qualifier: str,
        match_empty: bool = False,
    ):
        """
        Return True if the features share at least one common value for `qualifier`.
        `None` values in qualifiers are ignored.
        If `match_empty` is True then the features are connectable if they both
        do not have `qualifier` set.
        """
        f1q = [x for x in feature1.qualifiers.get(qualifier, []) if x is not None]
        f2q = [x for x in feature2.qualifiers.get(qualifier, []) if x is not None]
        if match_empty and not f1q and not f2q:
            return True
        else:
            isec = set(f1q).intersection(set(f2q))
            return len(isec) > 0

    @classmethod
    def _connectable_by_all_qualifiers(
        cls, feature1: SeqFeature, feature2: SeqFeature, by: str | Sequence[str]
    ):
        """
        Return True if `feature1` and `feature2` are connectable by all
        qualifiers in`by`.

        See also
        --------
        `_connectable_by_qualifier` method
        """
        if isinstance(by, str):
            by = (by,)
        ok = True
        for i, qualifier in enumerate(by):
            # If a top-level connection has been made allow missing qualifiers
            # to maintain the connection until an actual mismatch is
            # encountered
            match_empty = i > 0
            ok = cls._connectable_by_qualifier(
                feature1, feature2, qualifier, match_empty
            )
            if not ok:
                return False
        return ok

    @classmethod
    def _connectable_by_overlap(cls, feature1: SeqFeature, feature2: SeqFeature):
        """
        Return True if `feature1` fully overlaps `feature2`.

        See also
        --------
        gbtask.futils.covers
        """
        return futils.covers(feature1, feature2)

    def connect(self, by_qualifiers: str | Sequence[str] | None, by_overlap: bool):
        """
        Connect the feature graph by matching qualifier values and/or checking
        location overlaps.

        Parameters
        ----------

        by_qualifiers
            One or more feature qualifiers used to connect features. If None,
            features are not linked by qualifier.
        by_overlap
            If True, features are connected if they overlap. Overlap is only
            checked for features that were previously connected by qualifier.


        Returns
        -------

        The connected graph. Connected groups of features are connected
        components in the graph.
        """
        edges = list(self.graph.edges)
        if edges:
            # remove edges and start over
            self.graph.remove_edges_from(edges)
        if isinstance(by_qualifiers, str):
            by_qualifiers = (by_qualifiers,)

        for n1, n2 in itertools.combinations(self.graph.nodes, 2):
            if n1 == n2:
                continue
            _debug = f"Checking link {n1} -> {n2}"
            ok = True
            f1 = self.graph.nodes[n1]["feature"]
            f2 = self.graph.nodes[n2]["feature"]
            if by_qualifiers:
                ok = self._connectable_by_all_qualifiers(f1, f2, by_qualifiers)
            if ok and by_overlap:
                ok = ok and (
                    self._connectable_by_overlap(f1, f2)
                    or self._connectable_by_overlap(f2, f1)
                )
            if ok:
                _debug += f" - edge {(n1, n2)}"
                self.graph.add_edge(n1, n2)
            log.debug(_debug)
        return self.graph
