import logging

import pytest
from Bio.SeqFeature import CompoundLocation as cl
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation as sl
from gbtask.geneify.infer import ExonInferrer, UTRInferrer


class TestInferExons:
    def test_error_when_strands_differ(self):
        plus = SeqFeature(sl(1, 10, 1))
        minus = SeqFeature(sl(1, 10, -1))
        with pytest.raises(ValueError, match="Strands of all features must be equal"):
            ExonInferrer([plus, minus])

    @pytest.mark.parametrize(
        ("exon", "rna", "expect"),
        [
            (SeqFeature(sl(1, 10, 1)), SeqFeature(sl(1, 10, 1)), True),
            (SeqFeature(sl(1, 10, 1)), SeqFeature(sl(1, 20, 1)), True),
            (SeqFeature(sl(1, 10, -1)), SeqFeature(sl(1, 20, 1)), False),
            (SeqFeature(sl(1, 10, None)), SeqFeature(sl(1, 20, None)), True),
            (SeqFeature(sl(1, 10, 1)), SeqFeature(sl(15, 20, 1)), False),
            (SeqFeature(sl(20, 30, 1)), SeqFeature(sl(1, 10, 1)), False),
        ],
    )
    def test_is_compatible_with_simple_locations(self, exon, rna, expect, caplog):
        with caplog.at_level(logging.WARNING):
            assert ExonInferrer([]).is_compatible(exon, rna) is expect
            if expect is False:
                assert "is incompatible with RNA" in caplog.text

    @pytest.mark.parametrize(
        ("exon", "rna", "expect"),
        [
            (
                SeqFeature(sl(1, 10, 1)),
                SeqFeature(cl([sl(1, 10, 1), sl(20, 30, 1)])),
                True,
            ),
            (
                SeqFeature(sl(1, 10, 1)),
                SeqFeature(cl([sl(1, 20, 1), sl(30, 40, 1)])),
                True,
            ),
            (
                SeqFeature(sl(1, 10, -1)),
                SeqFeature(cl([sl(1, 20, 1), sl(1, 20, -1)])),
                False,
            ),
            (
                SeqFeature(sl(1, 10, 1)),
                SeqFeature(cl([sl(30, 40, 1), sl(50, 60, 1)])),
                False,
            ),
        ],
    )
    def test_is_compatible_with_compound_locations(self, exon, rna, expect, caplog):
        with caplog.at_level(logging.WARNING):
            assert ExonInferrer([]).is_compatible(exon, rna) is expect
            if expect is False:
                assert "is incompatible with RNA" in caplog.text

    def test__exons_from_rna(self):
        parts = ((1, 10), (20, 30), (40, 50))
        rna = SeqFeature(cl([sl(s, e) for s, e in parts]))
        expect = [SeqFeature(sl(s, e), type="exon") for s, e in parts]
        e = ExonInferrer([rna])
        assert e._exons_from_rna(rna) == expect

    def test__rna_from_exons(self):
        elocs = [sl(0, 10), sl(20, 30), sl(30, 40)]
        exons = [SeqFeature(l) for l in elocs]
        expect = SeqFeature(cl([sl(0, 10), sl(20, 40)]), type="someRNAtype")
        e = ExonInferrer(exons, rna="someRNAtype")
        assert e._rna_from_exons(exons) == [expect]

    def test_infer_prior_one_rna(self):
        rna = SeqFeature(cl([sl(0, 10), sl(20, 30)]), type="mRNA")
        elocs = [sl(0, 10), sl(20, 30)]
        exons = [SeqFeature(l, type="exon") for l in elocs]
        e = ExonInferrer([rna])
        assert e.infer() == exons

    def test_infer_prior_multiple_rnas(self):
        rna1 = SeqFeature(cl([sl(0, 10, 1), sl(20, 30, 1)]), type="mRNA")
        # shared exon with rna1
        rna2 = SeqFeature(cl([sl(0, 10, 1), sl(40, 50, 1)]), type="mRNA")
        elocs = [
            sl(0, 10, 1),  # rna1 and rna2
            sl(20, 30, 1),  # rna1
            sl(40, 50, 1),  # rna2
        ]
        exons = [SeqFeature(l, type="exon") for l in elocs]
        e = ExonInferrer([rna1, rna2])
        assert e.infer() == exons

    def test_infer_prior_exons(self):
        elocs = [sl(0, 10, 1), sl(20, 30, 1)]
        exons = [SeqFeature(l, type="exon") for l in elocs]
        rna = SeqFeature(cl([sl(0, 10, 1), sl(20, 30, 1)]), type="mRNA")
        e = ExonInferrer(exons)
        assert e.infer() == [rna]

    def test_infer_prior_rna_and_exons(self):
        rna1 = SeqFeature(cl([sl(0, 10, 1), sl(20, 30, 1), sl(40, 50, 1)]), type="mRNA")
        elocs_prior = [sl(0, 10, 1), sl(40, 50, 1)]
        exons_prior = [SeqFeature(l, type="exon") for l in elocs_prior]
        elocs_new = [sl(20, 30, 1)]
        exons_new = [SeqFeature(l, type="exon") for l in elocs_new]
        e = ExonInferrer(exons_prior + [rna1])
        assert e.infer() == exons_new


class TestInferUTRs:
    def test_error_when_strands_differ(self):
        plus = SeqFeature(sl(1, 10, 1))
        minus = SeqFeature(sl(1, 10, -1))
        with pytest.raises(ValueError, match="Strands of all features must be equal"):
            UTRInferrer([plus, minus])

    def test_1exon(self):
        exon = SeqFeature(sl(0, 100, 1), type="exon")
        cds = SeqFeature(sl(20, 80, 1), type="CDS")
        ui = UTRInferrer([exon, cds])
        expect = [(0, 20)], [(80, 100)]
        assert ui._utr_intervals() == expect

    def test_1exon_noutr(self):
        exon = SeqFeature(sl(0, 100, 1), type="exon")
        cds = SeqFeature(sl(0, 100, 1), type="CDS")
        ui = UTRInferrer([exon, cds])
        expect = ([], [])
        assert ui._utr_intervals() == expect

    def test_2exons(self):
        exon1 = SeqFeature(sl(0, 100, 1), type="exon")
        exon2 = SeqFeature(sl(200, 300, 1), type="exon")
        cds = SeqFeature(sl(20, 280, 1), type="CDS")
        ui = UTRInferrer([exon1, exon2, cds])
        expect = ([(0, 20)], [(280, 300)])
        assert ui._utr_intervals() == expect

    def test_2exons_noutr(self):
        exon1 = SeqFeature(sl(0, 100, 1), type="exon")
        exon2 = SeqFeature(sl(200, 300, 1), type="exon")
        cds = SeqFeature(sl(00, 300, 1), type="CDS")
        ui = UTRInferrer([exon1, exon2, cds])
        assert ui._utr_intervals() == ([], [])

    def test_2exons_noutr_cds_compoundloc(self):
        exon1 = SeqFeature(sl(0, 100, 1), type="exon")
        exon2 = SeqFeature(sl(200, 300, 1), type="exon")
        cds = SeqFeature(cl([sl(0, 100, 1), sl(200, 300, 1)]), type="CDS")
        ui = UTRInferrer([exon1, exon2, cds])
        assert ui._utr_intervals() == ([], [])

    def test_3exons_right_noncoding(self):
        exon1 = SeqFeature(sl(0, 100, 1), type="exon")
        exon2 = SeqFeature(sl(200, 300, 1), type="exon")
        exon3 = SeqFeature(sl(400, 500, 1), type="exon")
        cds = SeqFeature(sl(20, 280, 1), type="CDS")
        ui = UTRInferrer([exon1, exon2, exon3, cds])
        assert ui._utr_intervals() == ([(0, 20)], [(280, 300), (400, 500)])

    def test_3exons_left_noncoding(self):
        exon1 = SeqFeature(sl(0, 100, 1), type="exon")
        exon2 = SeqFeature(sl(200, 300, 1), type="exon")
        exon3 = SeqFeature(sl(400, 500, 1), type="exon")
        cds = SeqFeature(sl(220, 500, 1), type="CDS")
        ui = UTRInferrer([exon1, exon2, exon3, cds])
        assert ui._utr_intervals() == ([(0, 100), (200, 220)], [])

    def test_3exons_duplicated_unsorted(self):
        exon1 = SeqFeature(sl(0, 100, 1), type="exon")
        exon2 = SeqFeature(sl(200, 300, 1), type="exon")
        exon2_2 = SeqFeature(sl(200, 300, 1), type="exon")
        exon3 = SeqFeature(sl(400, 500, 1), type="exon")
        cds = SeqFeature(sl(220, 280, 1), type="CDS")
        ui = UTRInferrer([exon2_2, exon3, exon1, exon2, cds])
        assert ui._utr_intervals() == ([(0, 100), (200, 220)], [(280, 300), (400, 500)])

    def test_infer_plus(self):
        exon1 = SeqFeature(sl(0, 100, 1), type="exon")
        exon2 = SeqFeature(sl(200, 300, 1), type="exon")
        exon3 = SeqFeature(sl(400, 500, 1), type="exon")
        cds = SeqFeature(sl(220, 450, 1), type="CDS")
        ui = UTRInferrer([exon1, exon2, exon3, cds])
        utrs = ui.infer()
        expect = [
            SeqFeature(cl([sl(0, 100, 1), sl(200, 220, 1)]), type="5'UTR"),
            SeqFeature(sl(450, 500, 1), type="3'UTR"),
        ]
        assert utrs == expect

    def test_infer_minus(self):
        exon1 = SeqFeature(sl(0, 100, -1), type="exon")
        exon2 = SeqFeature(sl(200, 300, -1), type="exon")
        exon3 = SeqFeature(sl(400, 500, -1), type="exon")
        cds = SeqFeature(sl(220, 450, -1), type="CDS")
        ui = UTRInferrer([exon1, exon2, exon3, cds])
        utrs = ui.infer()
        expect = [
            SeqFeature(cl([sl(0, 100, -1), sl(200, 220, -1)]), type="3'UTR"),
            SeqFeature(sl(450, 500, -1), type="5'UTR"),
        ]
        assert utrs == expect

    def test_do_not_duplicate_existing(self):
        utr = SeqFeature(sl(1, 10, 1), type="5'UTR")
        exon = SeqFeature(sl(1, 30, 1), type="exon")
        cds = SeqFeature(sl(10, 20, 1), type="CDS")
        ui = UTRInferrer([utr, exon, cds])
        utrs = ui.infer()
        expect = [SeqFeature(sl(20, 30, 1), type="3'UTR")]
        assert utrs == expect
