import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from gbtask.togxf.GXFConverter import GFFConverter, GTFConverter, GXFConverter


@pytest.fixture
def record():
    r = SeqRecord(Seq(100 * "N"))
    r.name = "test"
    return r


@pytest.fixture
def gxfconverter(record):
    return GXFConverter(record)


@pytest.fixture
def gffconv(record):
    return GFFConverter(record)


@pytest.fixture
def gtfconv(record):
    return GTFConverter(record)


@pytest.fixture
def simpleCDS():
    f = SeqFeature(SimpleLocation(0, 100, 1), type="CDS")
    return f


@pytest.fixture
def compoundCDS():
    l1 = SimpleLocation(0, 100, 1)
    l2 = SimpleLocation(199, 300, 1)
    f = SeqFeature(CompoundLocation([l1, l2]), type="CDS")
    return f


@pytest.fixture
def simpleCDS_gff(simpleCDS):
    simpleCDS.qualifiers["ID"] = ["id1"]
    return simpleCDS


@pytest.fixture
def compoundCDS_gff(compoundCDS):
    compoundCDS.qualifiers["ID"] = ["id1"]
    return compoundCDS


@pytest.fixture
def simpleCDS_gtf(simpleCDS):
    simpleCDS.qualifiers["gene_id"] = ["id1"]
    return simpleCDS


@pytest.fixture
def compoundCDS_gtf(compoundCDS):
    compoundCDS.qualifiers["gene_id"] = ["id1"]
    return compoundCDS


class TestGXFConverter:
    def test__unsplittable(self, gxfconverter):
        assert gxfconverter._no_compound == {
            "gene",
            "mRNA",
            "snRNA",
            "snoRNA",
            "lncRNA",
            "tRNA",
            "ncRNA",
        }


class TestGFFConverter:
    def test_convert_type(self, gffconv):
        assert gffconv.convert_type("3'UTR") == "three_prime_utr"
        assert gffconv.convert_type("mRNA") == "mRNA"

    def test_fails_without_id(self, gffconv):
        f = SeqFeature(SimpleLocation(0, 10, 1), type="CDS")
        with pytest.raises(KeyError):
            gffconv.feature_to_lines(f)

    def test_fails_with_multiple_ids(self, gffconv):
        f = SeqFeature(
            SimpleLocation(0, 10, 1), type="CDS", qualifiers={"ID": ["id1", "id2"]}
        )
        with pytest.raises(KeyError):
            gffconv.feature_to_lines(f)

    def test_locations_simple(self, gffconv, simpleCDS_gff):
        result = gffconv._feature_locations(simpleCDS_gff)
        assert result == [("1", "100", "+")]

    def test_locations_compound(self, gffconv, compoundCDS_gff):
        result = gffconv._feature_locations(compoundCDS_gff)
        assert result == [("1", "100", "+"), ("200", "300", "+")]

    def test_locations_nocompound(self, gffconv, compoundCDS_gff):
        # change to a type that does not support compound locations in GFF
        compoundCDS_gff.type = "mRNA"
        result = gffconv._feature_locations(compoundCDS_gff)
        assert result == [("1", "300", "+")]

    def test_feature_to_lines_correct_number_of_lines(self, gffconv):
        l1 = SimpleLocation(0, 10, 1)
        l2 = SimpleLocation(20, 30, 1)
        f1 = SeqFeature(l1, qualifiers={"ID": "1"})
        assert len(gffconv.feature_to_lines(f1).split("\n")) == 1
        f2 = SeqFeature(CompoundLocation([l1, l2]), qualifiers={"ID": "1"})
        assert len(gffconv.feature_to_lines(f2).split("\n")) == 2
        f2 = SeqFeature(
            CompoundLocation([l1, l2]),
            qualifiers={"ID": "1", "Parent": ["1", "2"]},
        )
        assert len(gffconv.feature_to_lines(f2).split("\n")) == 2

    def test_feature_to_lines_simple(self, gffconv, simpleCDS_gff):
        expect = "\t".join(
            ["test", "gbtask_togff", "CDS", "1", "100", ".", "+", "0", "ID=id1"]
        )
        assert gffconv.feature_to_lines(simpleCDS_gff) == expect

    def test_feature_to_lines_compound(self, gffconv, compoundCDS_gff):
        line1 = "\t".join(
            ["test", "gbtask_togff", "CDS", "1", "100", ".", "+", "0", "ID=id1"]
        )
        line2 = "\t".join(
            # phase of 2!
            ["test", "gbtask_togff", "CDS", "200", "300", ".", "+", "2", "ID=id1"]
        )
        expect = "\n".join([line1, line2])
        assert gffconv.feature_to_lines(compoundCDS_gff) == expect


class TestGTFConverter:
    def test_convert_type(self, gtfconv):
        assert gtfconv.convert_type("3'UTR") == "UTR"
        assert gtfconv.convert_type("mRNA") == "mRNA"

    def test_fails_without_id(self, gtfconv):
        f = SeqFeature(SimpleLocation(0, 10, 1), type="CDS")
        with pytest.raises(KeyError):
            gtfconv.feature_to_lines(f)

    def test_fails_with_multiple_ids(self, gtfconv):
        f = SeqFeature(
            SimpleLocation(0, 10, 1), type="CDS", qualifiers={"gene_id": ["id1", "id2"]}
        )
        with pytest.raises(KeyError):
            gtfconv.feature_to_lines(f)

    def test_locations_simple(self, gtfconv, simpleCDS_gtf):
        result = gtfconv._feature_locations(simpleCDS_gtf)
        assert result == [("1", "100", "+")]

    def test_locations_compound(self, gtfconv, compoundCDS_gtf):
        result = gtfconv._feature_locations(compoundCDS_gtf)
        assert result == [("1", "100", "+"), ("200", "300", "+")]

    def test_locations_nocompound(self, gtfconv, compoundCDS_gtf):
        # change to a type that does not support compound locations in GFF
        compoundCDS_gtf.type = "mRNA"
        result = gtfconv._feature_locations(compoundCDS_gtf)
        assert result == [("1", "300", "+")]

    def test_feature_to_lines_correct_number_of_lines(self, gtfconv):
        l1 = SimpleLocation(0, 10, 1)
        l2 = SimpleLocation(20, 30, 1)
        f1 = SeqFeature(l1, qualifiers={"gene_id": "1"})
        assert len(gtfconv.feature_to_lines(f1).split("\n")) == 1
        f2 = SeqFeature(CompoundLocation([l1, l2]), qualifiers={"gene_id": "1"})
        assert len(gtfconv.feature_to_lines(f2).split("\n")) == 2
        f2 = SeqFeature(
            CompoundLocation([l1, l2]),
            qualifiers={"gene_id": "1", "transcript_id": ["1", "2"]},
        )
        assert len(gtfconv.feature_to_lines(f2).split("\n")) == 4

    def test_feature_to_lines_simple(self, gtfconv, simpleCDS_gtf):
        expect = "\t".join(
            ["test", "gbtask_togtf", "CDS", "1", "100", ".", "+", "0", 'gene_id "id1"']
        )
        assert gtfconv.feature_to_lines(simpleCDS_gtf) == expect

    def test_feature_to_lines_compound(self, gtfconv, compoundCDS_gtf):
        line1 = "\t".join(
            ["test", "gbtask_togtf", "CDS", "1", "100", ".", "+", "0", 'gene_id "id1"']
        )
        line2 = "\t".join(
            # phase of 2!
            [
                "test",
                "gbtask_togtf",
                "CDS",
                "200",
                "300",
                ".",
                "+",
                "2",
                'gene_id "id1"',
            ]
        )
        expect = "\n".join([line1, line2])
        assert gtfconv.feature_to_lines(compoundCDS_gtf) == expect

    def test_feature_to_lines_compound_multiple_transcript(
        self, gtfconv, compoundCDS_gtf
    ):
        txids = ["txid1", "txid2"]
        compoundCDS_gtf.qualifiers["transcript_id"] = txids
        expect = []
        for txid in txids:
            line1 = "\t".join(
                [
                    "test",
                    "gbtask_togtf",
                    "CDS",
                    "1",
                    "100",
                    ".",
                    "+",
                    "0",
                    f'gene_id "id1"; transcript_id "{txid}"',
                ]
            )
            line2 = "\t".join(
                # phase of 2!
                [
                    "test",
                    "gbtask_togtf",
                    "CDS",
                    "200",
                    "300",
                    ".",
                    "+",
                    "2",
                    f'gene_id "id1"; transcript_id "{txid}"',
                ]
            )
            expect.append(line1)
            expect.append(line2)
        expect = "\n".join(expect)
        assert gtfconv.feature_to_lines(compoundCDS_gtf) == expect

    def test_feature_to_lines_multiple_attribute_values(self, gtfconv, simpleCDS_gtf):
        simpleCDS_gtf.qualifiers["x"] = ["x1", "x2", "x3"]
        expect = []
        expect = "\t".join(
            [
                "test",
                "gbtask_togtf",
                "CDS",
                "1",
                "100",
                ".",
                "+",
                "0",
                'gene_id "id1"; x "x1"; x "x2"; x "x3"',
            ]
        )
        assert gtfconv.feature_to_lines(simpleCDS_gtf) == expect
