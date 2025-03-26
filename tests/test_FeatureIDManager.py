from itertools import permutations

import pytest
from gbtask.utils.FeatureIDManager import FeatureIDManager

from .helpers import feature


@pytest.mark.parametrize(
    ("features", "prefix"),
    [
        ([feature("exon", expect="exon.01")], ""),
        ([feature("exon", expect="exon.01"), feature("exon", expect="exon.02")], ""),
        (
            [
                feature("exon", expect="prefix:exon.01"),
                feature("exon", expect="prefix:exon.02"),
            ],
            "prefix",
        ),
        (
            [
                feature("exon", expect="prefix:exon.01"),
                feature("exon", expect="prefix:exon.02"),
                feature("mRNA", expect="prefix:mRNA.01"),
            ],
            "prefix",
        ),
    ],
)
def test_set_feature_id(features, prefix):
    mgr = FeatureIDManager(features, prefix)
    mgr.set_feature_ids()
    for feature in mgr.features:
        expect = feature.qualifiers["expect"]
        assert feature.id == expect
        assert feature.qualifiers["ID"] == [expect]


def test_set_id_for_equal_features():
    f1 = feature("gene")
    f2 = feature("gene")
    mgr = FeatureIDManager([f1, f2])
    mgr.set_feature_ids()
    assert f1.id == "gene.01"
    assert f2.id == "gene.02"


def test_set_id_for_duplicate_features():
    f1 = feature("gene")
    f2 = f1
    mgr = FeatureIDManager([f1, f2])
    mgr.set_feature_ids()
    assert f1 is f2
    assert f1.id == "gene.01"
    assert f1.id == f2.id


@pytest.mark.parametrize("order", permutations([0, 1, 2]))
def test_set_feature_ids_sorts(order):
    exon1 = feature("exon", 0, 10)
    exon2 = feature("exon", 10, 20)
    exon3 = feature("exon", 20, 30)
    l = [exon1, exon2, exon3]
    # reorder
    l = [l[i] for i in order]
    mgr = FeatureIDManager(l)
    mgr.set_feature_ids()
    assert exon1.id == "exon.01"
    assert exon2.id == "exon.02"
    assert exon3.id == "exon.03"
