# import logging
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation as cl
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation as sl
from Bio.SeqRecord import SeqRecord

from gbtask.splice.splice import _truncate, splice


@pytest.fixture
def record():
    return SeqRecord(Seq("AATTGGCCAT"))


def test_truncate_strand_mismatch():
    covered = SeqFeature(sl(0, 10, 1))
    to = SeqFeature(sl(0, 10, -1))
    assert _truncate(covered, to, stranded=True) is None


def test_truncate_strand_match():
    covered = SeqFeature(sl(0, 10, 1))
    to = SeqFeature(sl(3, 7, 1))
    assert _truncate(covered, to, stranded=True) == to


def test_truncate_spanning():
    covered = SeqFeature(sl(0, 10))
    to = SeqFeature(sl(0, 10))
    assert _truncate(covered, to) == covered


def test_truncate_within():
    covered = SeqFeature(sl(1, 9))
    to = SeqFeature(sl(0, 10))
    assert _truncate(covered, to) == covered


def test_truncate_exceeds_left():
    covered = SeqFeature(sl(0, 10))
    to = SeqFeature(sl(1, 10))
    assert _truncate(covered, to) == SeqFeature(sl(1, 10))


def test_truncate_exceeds_right():
    covered = SeqFeature(sl(0, 11))
    to = SeqFeature(sl(0, 10))
    assert _truncate(covered, to) == SeqFeature(sl(0, 10))


def test_truncate_exceeds_both():
    covered = SeqFeature(sl(0, 11))
    to = SeqFeature(sl(1, 10))
    assert _truncate(covered, to) == SeqFeature(sl(1, 10))


def test_splice_0exons(record):
    spliced = splice(record, "exon")
    assert spliced.seq == ""


def test_splice_1exon(record):
    exon = SeqFeature(sl(0, 2), type="exon")
    record.features.append(exon)
    spliced = splice(record, "exon")
    assert spliced.seq == "AA"
    assert spliced.features == [exon]


def test_splice_2exons(record):
    exon1 = SeqFeature(sl(0, 2), type="exon")
    exon2 = SeqFeature(sl(4, 6), type="exon")
    record.features = [exon1, exon2]
    spliced = splice(record, "exon")
    exon2_aftersplice = SeqFeature(sl(2, 4), type="exon")
    assert spliced.seq == "AAGG"
    assert spliced.features == [exon1, exon2_aftersplice]


def test_splice_2exons_adjacent(record):
    exon1 = SeqFeature(sl(0, 2), type="exon")
    exon2 = SeqFeature(sl(2, 4), type="exon")
    record.features = [exon1, exon2]
    spliced = splice(record, "exon")
    assert spliced.seq == "AATT"
    assert spliced.features == [exon1, exon2]


def test_splice_2exons_overlapping(record):
    exon1 = SeqFeature(sl(0, 2), type="exon")
    exon2 = SeqFeature(sl(1, 3), type="exon")
    record.features = [exon1, exon2]
    spliced = splice(record, "exon")
    assert spliced.seq == "AAT"
    assert spliced.features == [exon1, exon2]


def test_splice_keeps_covered_simplefeatures(record):
    exon = SeqFeature(sl(0, 5), type="exon")
    other = SeqFeature(sl(3, 5), type="other")
    record.features = [exon, other]
    spliced = splice(record, "exon")
    other_aftersplice = SeqFeature(sl(3, 5), type="other")
    assert spliced.features == [exon, other_aftersplice]


def test_splice_keeps_covered_compoundfeatures(record):
    exon1 = SeqFeature(sl(0, 5), type="exon")
    exon2 = SeqFeature(sl(7, 9), type="exon")
    other = SeqFeature(cl([sl(3, 5), sl(7, 9)]), type="other")
    record.features = [exon1, exon2, other]
    spliced = splice(record, "exon")
    exon2_aftersplice = SeqFeature(sl(5, 7), type="exon")
    other_aftersplice1 = SeqFeature(sl(3, 5), type="other")
    other_aftersplice2 = SeqFeature(sl(5, 7), type="other")
    assert spliced.features == [
        exon1,
        other_aftersplice1,
        exon2_aftersplice,
        other_aftersplice2,
    ]


def test_splice_truncates_and_splits_features(record):
    exon1 = SeqFeature(sl(0, 5), type="exon")
    exon2 = SeqFeature(sl(7, 10), type="exon")
    other = SeqFeature(sl(0, 12), type="other")
    record.features = [exon1, exon2, other]
    spliced = splice(record, "exon")
    exon2_aftersplice = SeqFeature(sl(5, 8), type="exon")
    other_aftersplice1 = SeqFeature(sl(0, 5), type="other")
    other_aftersplice2 = SeqFeature(sl(5, 8), type="other")
    assert spliced.features == [
        exon1,
        other_aftersplice1,
        exon2_aftersplice,
        other_aftersplice2,
    ]


def test_splice_truncates_and_splits_features_with_partial_overlap(record):
    exon1 = SeqFeature(sl(0, 5), type="exon")
    exon2 = SeqFeature(sl(7, 10), type="exon")
    other = SeqFeature(sl(2, 8), type="other")
    record.features = [exon1, exon2, other]
    spliced = splice(record, "exon")
    exon2_aftersplice = SeqFeature(sl(5, 8), type="exon")
    other_aftersplice1 = SeqFeature(sl(2, 5), type="other")
    other_aftersplice2 = SeqFeature(sl(5, 6), type="other")
    assert spliced.features == [
        exon1,
        other_aftersplice1,
        exon2_aftersplice,
        other_aftersplice2,
    ]
