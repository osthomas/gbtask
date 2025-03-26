import logging
from typing import Hashable, Sequence, Set

from Bio.SeqRecord import SeqRecord

from gbtask.geneify.infer import ExonInferrer, Inferrer, UTRInferrer
from gbtask.utils.FeatureGroups import FeatureGroups
from gbtask.utils.FeatureHierarchy import FeatureHierarchy
from gbtask.utils.futils import split_by_strand

log = logging.getLogger(__name__)


def geneify(
    record: SeqRecord,
    group_by: Sequence[str],
    byoverlap: bool,
    infer_exons: bool,
    infer_utrs: bool,
):
    features = [f for f in record.features]
    fg = FeatureGroups(features)
    fg.connect(by_qualifiers=group_by, by_overlap=byoverlap)
    # Perform feature inference per group
    if infer_exons or infer_utrs:
        if infer_exons:
            # TODO: inferred RNAs are always 'mRNA', allow detection
            # of other RNA types per group
            fg.per_group(_infer_exons)
        if infer_utrs:
            fg.per_group(_infer_utrs)
    fh = FeatureHierarchy(fg.graph, add_root=True)
    fh.to_features()
    record.features = fh.features
    return record


def _infer_exons(featuregroups: FeatureGroups, nodes: Set[Hashable]):
    _infer_per_strand(featuregroups, nodes, ExonInferrer)


def _infer_utrs(featuregroups: FeatureGroups, nodes: Set[Hashable]):
    _infer_per_strand(featuregroups, nodes, UTRInferrer)


def _infer_per_strand(
    featuregroups: FeatureGroups, nodes: Set[Hashable], inferrer: type[Inferrer]
):
    features = [featuregroups.graph.nodes[node]["feature"] for node in nodes]
    for strand_features in split_by_strand(features).values():
        new_features = inferrer(strand_features).infer()
        for new_feature in new_features:
            node_id = featuregroups.add_feature_node(new_feature)
            # Connect to one node of the group so the new feature becomes part
            # of the group. The exact connection does not matter.
            featuregroups.graph.add_edge(next(iter(nodes)), node_id)
