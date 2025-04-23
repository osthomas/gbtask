import pytest

from gbtask.utils.TypeTree import TypeTree


@pytest.fixture
def tree():
    return TypeTree()


def test_root_type(tree):
    assert tree.root_type == "gene"


def test_alias(tree):
    assert tree.alias("3'UTR", "gff") == "three_prime_UTR"
    assert tree.alias("3'UTR", "nonexistent") == "3'UTR"


def test_relabel_gff(tree):
    want_orig = ("5'UTR", "3'UTR")
    want_new = ("five_prime_UTR", "three_prime_UTR")
    for wo, wn in zip(want_orig, want_new):
        assert wo in tree._tree.nodes
        assert wn not in tree._tree.nodes
    g2 = tree.relabel("gff")
    for wo, wn in zip(want_orig, want_new):
        assert wo not in g2._tree.nodes
        assert wn in g2._tree.nodes
    # original should be unchanged
    for wo, wn in zip(want_orig, want_new):
        assert wo in tree._tree.nodes
        assert wn not in tree._tree.nodes


def test_relabel_prevent_duplicates(tree):
    tree._tree.nodes["exon"]["alias"] = {"gff": "three_prime_UTR"}
    with pytest.raises(NotImplementedError, match="not unique"):
        tree.relabel("gff")


def test_parent_types(tree):
    assert not list(tree.parent_types("gene"))
    assert list(tree.parent_types("mRNA")) == ["gene"]
    assert set(tree.parent_types("exon")) == {
        "mRNA",
        "tRNA",
        "snoRNA",
        "snRNA",
        "ncRNA",
        "lncRNA",
    }


def test_parent_type_of_nonexistent(tree):
    assert tree.parent_types("atcgatcga") == tuple()


def test_gtf_levels(tree):
    assert tree.gtf_levels == ["gene_id", "transcript_id"]


def test_gtf_level_types(tree):
    assert tree.gtf_level_types("gene_id") == ["gene"]
    assert tree.gtf_level_types("transcript_id") == [
        "mRNA",
        "tRNA",
        "snRNA",
        "snoRNA",
        "lncRNA",
        "ncRNA",
    ]
