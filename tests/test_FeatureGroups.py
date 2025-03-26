import networkx as nx
import pytest
from gbtask.utils.FeatureGroups import FeatureGroups

from .helpers import feature
from .helpers import features as features_


@pytest.fixture
def features():
    return features_()


def test_init():
    fg = FeatureGroups([feature("gene"), feature("mRNA"), feature("CDS")])
    assert len(fg.graph) == 3
    assert len(fg.graph.edges) == 0


def test_init_error_if_id_not_unique():
    f1 = feature()
    f1.id = "X"
    f2 = feature()
    f2.id = "X"
    with pytest.raises(ValueError):
        FeatureGroups([f1, f2], node_label_from_id=True)


@pytest.mark.parametrize(
    ("qualifiers1", "qualifiers2", "expect"),
    [
        ({"gene": []}, {"gene": []}, False),
        ({"gene": []}, {"gene": ["A"]}, False),
        ({"gene": ["A"]}, {"gene": ["A"]}, True),
        ({"gene": ["A2"]}, {"gene": ["A"]}, False),
        ({"gene": ["A"]}, {"gene": ["B"]}, False),
        ({"gene": ["A", "B"]}, {"gene": ["B"]}, True),
        ({"gene": ["A", "B"]}, {"gene": ["B", "C"]}, True),
        ({"gene": ["A", "B"]}, {"gene": ["B", "C"], "transcript": ["A"]}, True),
    ],
)
def test_connectable_by_qualifier(qualifiers1, qualifiers2, expect):
    f1 = feature(**qualifiers1)
    f2 = feature(**qualifiers2)
    assert FeatureGroups._connectable_by_qualifier(f1, f2, "gene") == expect


@pytest.mark.parametrize(
    ("qualifiers1", "qualifiers2", "expect"),
    [
        ({"gene": []}, {"gene": []}, True),
        ({"gene": []}, {"gene": ["A"]}, False),
        ({"gene": ["A"]}, {"gene": ["A"]}, True),
        ({"gene": ["A2"]}, {"gene": ["A"]}, False),
        ({"gene": ["A"]}, {"gene": ["B"]}, False),
        ({"gene": ["A", "B"]}, {"gene": ["B"]}, True),
        ({"gene": ["A", "B"]}, {"gene": ["B", "C"]}, True),
        ({"gene": ["A", "B"]}, {"gene": ["B", "C"], "transcript": ["A"]}, True),
    ],
)
def test_connectable_by_qualifier_match_empty(qualifiers1, qualifiers2, expect):
    f1 = feature(**qualifiers1)
    f2 = feature(**qualifiers2)
    assert (
        FeatureGroups._connectable_by_qualifier(f1, f2, "gene", match_empty=True)
        == expect
    )


@pytest.mark.parametrize(
    ("qualifiers1", "qualifiers2", "by", "expect"),
    [
        ({"gene": []}, {"gene": []}, "gene", False),
        ({"gene": []}, {"gene": ["A"]}, "gene", False),
        ({"gene": ["A"]}, {"gene": ["A"]}, "gene", True),
        ({"gene": ["A2"]}, {"gene": ["A"]}, "gene", False),
        ({"gene": ["A"]}, {"gene": ["B"]}, "gene", False),
        ({"gene": ["A", "B"]}, {"gene": ["B"]}, "gene", True),
        ({"gene": ["A", "B"]}, {"gene": ["B", "C"]}, "gene", True),
        ({"gene": ["A", "B"]}, {"gene": ["B", "C"], "tx": ["A"]}, "gene", True),
        (
            {"gene": ["A", "B"]},
            {"gene": ["B", "C"], "tx": ["A"]},
            "tx",
            False,
        ),
        (
            {"gene": ["A", "B"], "tx": ["A"]},
            {"gene": ["B", "C"], "tx": ["A"]},
            "tx",
            True,
        ),
        (
            {"gene": ["A"]},
            {"gene": ["A"]},
            ("gene", "tx"),
            True,
        ),
        (
            {"gene": ["A"], "tx": []},
            {"gene": ["A"]},
            ("gene", "tx"),
            True,
        ),
        (
            {"gene": ["A"], "tx": ["A"]},
            {"gene": ["A"], "tx": ["A"]},
            ("gene", "tx"),
            True,
        ),
        (
            {"gene": ["A", "B"], "tx": ["A"]},
            {"gene": ["A"], "tx": ["A"]},
            ("gene", "tx"),
            True,
        ),
        (
            {"gene": ["A", "B"], "tx": ["A"], "prod": ["X"]},
            {"gene": ["A"], "tx": ["A"]},
            ("gene", "tx"),
            True,
        ),
        (
            {"gene": ["A", "B"], "tx": ["A"], "prod": ["X"]},
            {"gene": ["A"], "tx": ["A"]},
            ("gene", "tx", "prod"),
            False,
        ),
        (
            {"gene": ["A", "B"], "tx": [], "prod": ["X"]},
            {"gene": ["A"], "prod": ["X"]},
            ("gene", "tx", "prod"),
            True,
        ),
    ],
)
def test_connectable_by_all_qualifiers(qualifiers1, qualifiers2, by, expect):
    f1 = feature(**qualifiers1)
    f2 = feature(**qualifiers2)
    assert FeatureGroups._connectable_by_all_qualifiers(f1, f2, by) == expect


@pytest.mark.parametrize(
    ("feature1", "feature2", "expect"),
    [
        (feature(start=0, end=10), feature(start=0, end=10), True),
        (feature(start=1, end=10), feature(start=0, end=10), False),
    ],
)
def test_connectable_by_overlap(feature1, feature2, expect):
    assert FeatureGroups._connectable_by_overlap(feature1, feature2) == expect


def test_connect_qualifiers_and_overlap(features):
    fg = FeatureGroups(features.values(), node_label_from_id=True)
    fg.connect("gene", True)
    assert nx.number_connected_components(fg.graph) == 5
    for gene in ("A", "B", "A2"):
        gene_nodes = set([k for k in features if k.endswith(gene)])
        assert nx.node_connected_component(fg.graph, "gene" + gene) == gene_nodes
    assert nx.node_connected_component(fg.graph, "dangling1") == {"dangling1"}
    assert nx.node_connected_component(fg.graph, "dangling2") == {"dangling2"}


def test_connect_qualifiers_only(features):
    fg = FeatureGroups(features.values(), node_label_from_id=True)
    fg.connect("gene", False)
    fg.connect("gene", False)
    assert nx.number_connected_components(fg.graph) == 4
    for gene in ("B",):
        gene_nodes = set([k for k in features if k.endswith(gene)])
        assert nx.node_connected_component(fg.graph, "gene" + gene) == gene_nodes
    assert nx.node_connected_component(fg.graph, "geneA") == {
        "geneA",
        "mrnaA",
        "exon1A",
        "exon2A",
        "cdsA",
        "geneA2",
        "mrnaA2",
        "exon1A2",
        "exon2A2",
        "cdsA2",
    }
    assert nx.node_connected_component(fg.graph, "dangling1") == {"dangling1"}
    assert nx.node_connected_component(fg.graph, "dangling2") == {"dangling2"}


def test_connect_overlap_only(features):
    fg = FeatureGroups(features.values(), node_label_from_id=True)
    fg.connect(None, True)
    assert nx.number_connected_components(fg.graph) == 3
    assert nx.node_connected_component(fg.graph, "geneA") == {
        "dangling1",
        "geneA",
        "mrnaA",
        "exon1A",
        "exon2A",
        "cdsA",
        "geneB",
        "mrnaB",
        "exon1B",
        "exon2B",
        "cdsB",
    }
    assert nx.node_connected_component(fg.graph, "geneA2") == {
        "geneA2",
        "mrnaA2",
        "exon1A2",
        "exon2A2",
        "cdsA2",
    }
    assert nx.node_connected_component(fg.graph, "dangling2") == {"dangling2"}


def test_connect_repeatedly(features):
    fg = FeatureGroups(features.values(), node_label_from_id=True)
    fg.connect("gene", True)
    # running again recreates internal graph
    fg.connect("gene", True)
    in_features = list(features.values())
    # Features should refer to the original objects after re-creation of the graph
    graph_features = []
    for n in fg.graph.nodes:
        graph_features.append(fg.graph.nodes[n]["feature"])
    assert len(in_features) == len(graph_features)
    for gf in graph_features:
        assert gf in in_features


# def test_validation_fails_with_0_roots():
#     ft = FeatureTree([mrna(0, 10)], TypeTree)
#     with pytest.raises(FeatureTreeError, match="but got 0"):
#         ft._validate()


# def test_validation_fails_with_2_roots():
#     ft = FeatureTree(gene("A") + gene("B"), TypeTree)
#     with pytest.raises(FeatureTreeError, match="but got 2"):
#         ft._validate()


# def test_get_correct_root_and_leaf_features():
#     gene = feature("gene")
#     mrna = feature("mRNA")
#     cds = feature("CDS")
#     exon = feature("exon")
#     random = feature("random")
#     ft = FeatureTree([gene, mrna, cds, exon, random], TypeTree)
#     assert ft.root_features == [gene]
#     assert ft.leaf_features == [cds, exon]


# def test_get_parent_features(caplog):
#     gene = feature("gene")
#     mrna = feature("mRNA")
#     cds = feature("CDS")
#     exon = feature("exon")
#     random = feature("random")
#     ft = FeatureTree([gene, mrna, cds, exon, random], TypeTree)
#     assert ft._get_parents(gene) == []
#     assert ft._get_parents(mrna) == [gene]
#     assert ft._get_parents(cds) == [mrna]
#     assert ft._get_parents(exon) == [mrna]
#     with caplog.at_level(logging.INFO):
#         assert ft._get_parents(random) == []
#         assert "random is not in the feature type hierarchy" in caplog.text
