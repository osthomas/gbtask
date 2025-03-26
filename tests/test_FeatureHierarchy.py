import networkx as nx
import pytest
from Bio.SeqFeature import SeqFeature
from gbtask.utils.FeatureGroups import FeatureGroups
from gbtask.utils.FeatureHierarchy import FeatureHierarchy

from .helpers import feature
from .helpers import features as features_


@pytest.fixture
def empty():
    return FeatureHierarchy(nx.Graph())


def prepare_featuregraph(features: list[SeqFeature]):
    """
    Helper to instantiate FeatureGroups with proper labels to make debugging
    easier.
    """
    fg = FeatureGroups(features, node_label_from_id=True)
    fg.connect("gene", True)
    return fg


def prepare_featurehierarchy(features: list[SeqFeature]):
    fg = prepare_featuregraph(features)
    fh = FeatureHierarchy(fg.graph)
    return fh


@pytest.fixture
def connected():
    features = list(features_().values())
    return prepare_featuregraph(features)


@pytest.fixture
def hierarchy():
    features = list(features_().values())
    fg = prepare_featuregraph(features)
    return FeatureHierarchy(fg.graph)


def hierarchies_equal(
    result: nx.DiGraph, expect: dict[str, dict], visual: bool = False
):
    # generate `expect` dict via `nx.to_dict_of_dicts(result)` after visual assertion
    expect_graph = nx.from_dict_of_dicts(
        expect,
        create_using=nx.DiGraph,
    )
    # steal features from expected result
    for node in expect_graph.nodes:
        try:
            expect_graph.nodes[node]["feature"] = result.nodes[node]["feature"]
        except KeyError:
            # will be caught my graph_equal
            pass
    if visual:
        print(nx.to_dict_of_dicts(result))
        print(nx.to_dict_of_dicts(expect_graph))
        # Visual assertion:
        from matplotlib import pyplot as plt

        pos = nx.planar_layout(result)
        nx.draw_networkx(result, pos=pos)
        plt.show()

    ok = nx.utils.graphs_equal(result, expect_graph)
    return ok


def test_init(empty):
    pass


def test_node_groups(hierarchy):
    assert hierarchy.node_groups == [
        {"geneA", "mrnaA", "exon1A", "exon2A", "cdsA"},
        {"geneB", "mrnaB", "exon1B", "exon2B", "cdsB"},
        {"geneA2", "mrnaA2", "exon1A2", "exon2A2", "cdsA2"},
        {"dangling1"},
        {"dangling2"},
    ]


def test_feature_groups():
    fh = FeatureHierarchy(nx.DiGraph())
    f1 = feature()
    f2 = feature()
    f3 = feature()
    fh.graph.add_node(1, feature=f1)
    fh.graph.add_node(2, feature=f2)
    fh.graph.add_node(3, feature=f3)
    fh.graph.add_edge(2, 3)
    assert fh.feature_groups == [
        [f1],
        [f2, f3],
    ]


def test_is_parent(empty):
    assert empty._is_parent(feature("gene"), feature("mRNA")) is True
    assert empty._is_parent(feature("mRNA"), feature("exon")) is True
    assert empty._is_parent(feature("exon"), feature("mRNA")) is False


def test_hierarchy_keeps_features(hierarchy):
    for node in hierarchy.graph.nodes:
        assert isinstance(hierarchy.graph.nodes[node]["feature"], SeqFeature)


def test_hierarchize_reproducible(hierarchy):
    # Make copys via `to_dict_of_dicts` because graphs are mutable
    prev = nx.to_dict_of_dicts(hierarchy.graph)
    for _ in range(3):
        hierarchy._hierarchize()
        this = nx.to_dict_of_dicts(hierarchy.graph)
        assert prev == this


def test_hierarchize(hierarchy):
    result = hierarchy.graph
    expect = {
        "mrnaA": {"exon2A": {}, "cdsA": {}, "exon1A": {}},
        "exon2A": {},
        "cdsA": {},
        "exon1A": {},
        "geneA": {"mrnaA": {}},
        "cdsB": {},
        "geneB": {"mrnaB": {}},
        "mrnaB": {"cdsB": {}, "exon2B": {}, "exon1B": {}},
        "exon2B": {},
        "exon1B": {},
        "cdsA2": {},
        "mrnaA2": {"cdsA2": {}, "exon2A2": {}, "exon1A2": {}},
        "geneA2": {"mrnaA2": {}},
        "exon2A2": {},
        "exon1A2": {},
        "dangling1": {},
        "dangling2": {},
    }
    assert hierarchies_equal(result, expect)


def test_hierarchize_multiple_roots():
    graph = nx.Graph()
    graph.add_node("geneA", feature=feature("gene"))
    graph.add_node("geneB", feature=feature("gene"))
    graph.add_node("CDS", feature=feature("CDS"))
    graph.add_edge("geneA", "CDS")
    graph.add_edge("geneB", "CDS")
    fh = FeatureHierarchy(graph)
    fh._hierarchize()
    expect = {"geneA": {"CDS": {}}, "geneB": {"CDS": {}}, "CDS": {}}
    assert hierarchies_equal(fh.graph, expect)


def test_hierarchize_skip_missing_intermediate():
    graph = nx.Graph()
    graph.add_node("geneX", feature=feature("gene"))
    graph.add_node("CDSX", feature=feature("CDS"))
    graph.add_edge("geneX", "CDSX")
    fh = FeatureHierarchy(graph)
    fh._hierarchize()
    expect = {"geneX": {"CDSX": {}}, "CDSX": {}}
    assert hierarchies_equal(fh.graph, expect)


def test_hierarchize_link_unknown_to_root():
    graph = nx.Graph()
    graph.add_node("geneX", feature=feature("gene"))
    graph.add_node("random", feature=feature("random"))
    graph.add_edge("geneX", "random")
    fh = FeatureHierarchy(graph)
    fh._hierarchize()
    expect = {"geneX": {"random": {}}, "random": {}}
    assert hierarchies_equal(fh.graph, expect)


def test_add_root_node_to():
    f = feature("mRNA")
    fh = prepare_featurehierarchy([f])
    node_id = fh._add_root_node_to(set(fh.graph))
    assert fh.graph.nodes[node_id]["feature"].type == "gene"


def test_add_root_node_to_dont_add_for_root_type():
    f = feature("gene")
    fh = prepare_featurehierarchy([f])
    assert fh._add_root_node_to(set(fh.graph)) is None


def test_add_root_node_to_dont_add_if_root_type_in_group():
    graph = nx.DiGraph()
    f1 = feature("gene")
    f2 = feature("mRNA")
    graph.add_node("geneA", feature=f1)
    graph.add_node("mrnaA", feature=f2)
    graph.add_edge("geneA", "mrnaA")
    fh = FeatureHierarchy(graph)
    fh.graph = graph
    assert fh._add_root_node_to(set(fh.graph)) is None


def test_add_root_node_to_dont_add_for_stray_type():
    f = feature("notinhierarchy")
    fh = prepare_featurehierarchy([f])
    assert fh._add_root_node_to(set(fh.graph)) is None


def test_add_root_node_to_stray_features_dont_enlarge_newroot():
    f1 = feature("notinhierarchy", 100, 200, id="0")
    f2 = feature("mRNA", 0, 50, id="1")
    fh = prepare_featurehierarchy([f1, f2])
    node_id = fh._add_root_node_to(set(fh.graph))
    assert len(fh.graph) == 3
    newroot = fh.graph.nodes[node_id]["feature"]
    assert newroot.type == "gene"
    assert newroot.location.start == 0
    assert newroot.location.end == 50


def test_hierarchize_with_add_root():
    f1 = feature("mRNA", gene=["geneA"], id="0")
    f2 = feature("mRNA", gene=["geneB"], id="1")
    f3 = feature("misc_feature", id="2")
    fh = prepare_featurehierarchy([f1, f2, f3])
    fh._hierarchize(add_root=True)
    assert len(fh.graph) == 5
    for feature_group in fh.feature_groups:
        if f1 in feature_group:
            assert len([f for f in feature_group if f.type == "gene"]) == 1
        if f2 in feature_group:
            assert len([f for f in feature_group if f.type == "gene"]) == 1
        if f3 in feature_group:
            assert len([f for f in feature_group if f.type == "gene"]) == 0


def test_set_feature_ids_method(hierarchy):
    for feature in hierarchy.features:
        feature.__before = feature.id
    hierarchy.set_feature_ids()
    before2after = {
        "geneA": "gene.01",
        "mrnaA": "gene.01:mRNA.01",
        "exon1A": "gene.01:exon.01",
        "exon2A": "gene.01:exon.02",
        "cdsA": "gene.01:CDS.01",
        "geneB": "gene.02",
        "mrnaB": "gene.02:mRNA.01",
        "exon1B": "gene.02:exon.01",
        "exon2B": "gene.02:exon.02",
        "cdsB": "gene.02:CDS.01",
        "geneA2": "gene.03",
        "mrnaA2": "gene.03:mRNA.01",
        "exon1A2": "gene.03:exon.01",
        "exon2A2": "gene.03:exon.02",
        "cdsA2": "gene.03:CDS.01",
        "dangling1": "misc_feature.01",
        "dangling2": "regulatory.01",
    }
    assert len(before2after) == len(hierarchy.features)
    for feature in hierarchy.features:
        assert before2after[feature.__before] == feature.id


def test_feature_ids_method_feature_in_multiple_groups():
    features = [
        feature("gene", gene=["geneA"], expect=["gene.01"]),
        feature("gene", gene=["geneB"], expect=["gene.02"]),
        feature("exon", gene=["geneA", "geneB"], expect=["gene.01_gene.02:exon.01"]),
    ]
    fg = FeatureGroups(features)
    graph = fg.connect("gene", True)
    fh = FeatureHierarchy(graph)
    fh.set_feature_ids()
    for f in fh.features:
        assert f.id == f.qualifiers["expect"][0]


def test_feature_ids_method_feature_with_multiple_roots():
    features = [
        feature("gene", gene=["geneA"], expect=["gene.01"]),
        feature("gene", gene=["geneA"], expect=["gene.02"]),
        feature("exon", gene=["geneA", "geneB"], expect=["gene.01_gene.02:exon.01"]),
    ]
    fg = FeatureGroups(features)
    graph = fg.connect("gene", True)
    fh = FeatureHierarchy(graph)
    fh.set_feature_ids()
    for f in fh.features:
        assert f.id == f.qualifiers["expect"][0]


def test_set_gff_qualifiers():
    features = [
        feature("gene", id="geneA", gene=["geneA"], expect=[]),
        feature("mRNA", id="mrnaA1", gene=["geneA"], expect=["geneA"]),
        feature("mRNA", id="mrnaA2", gene=["geneA"], expect=["geneA"]),
        feature("exon", id="exonA1", gene=["geneA"], expect=["mrnaA1", "mrnaA2"]),
    ]
    hierarchy = prepare_featurehierarchy(features)
    hierarchy._set_gff_qualifiers()
    for node in hierarchy.graph:
        f = hierarchy.graph.nodes[node]["feature"]
        assert f.qualifiers["Parent"] == f.qualifiers["expect"]


def test_set_gtf_qualifiers():
    features = [
        feature(
            "gene",
            id="geneA",
            gene=["geneA"],
            expect_gene=["geneA"],
            expect_transcript=[],
        ),
        feature(
            "mRNA",
            id="mrnaA1",
            gene=["geneA"],
            expect_gene=["geneA"],
            expect_transcript=["mrnaA1"],
        ),
        feature(
            "mRNA",
            id="mrnaA2",
            gene=["geneA"],
            expect_gene=["geneA"],
            expect_transcript=["mrnaA2"],
        ),
        feature(
            "exon",
            id="exonA1",
            gene=["geneA"],
            expect_gene=["geneA"],
            expect_transcript=["mrnaA1", "mrnaA2"],
        ),
    ]
    hierarchy = prepare_featurehierarchy(features)
    hierarchy._set_gtf_qualifiers()
    for node in hierarchy.graph:
        f = hierarchy.graph.nodes[node]["feature"]
        assert f.qualifiers["gene_id"] == f.qualifiers["expect_gene"]
        if f.type == "gene":
            # genes have no transcript_id
            assert "transcript_id" not in f.qualifiers
        else:
            assert f.qualifiers["transcript_id"] == f.qualifiers["expect_transcript"]
