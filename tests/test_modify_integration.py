import subprocess
from importlib.resources import files

import pytest

from tests import integration


def case_io(prefix: str):
    """
    Get input and output files for a test case, assuming a common prefix.
    """
    import tests

    inf = files(tests).joinpath("resources", "integration_modify", prefix + "_input.gb")
    outf = files(tests).joinpath(
        "resources", "integration_modify", prefix + "_output.gb"
    )
    return inf, outf


def modify_run(args, infile):
    """
    Run the modify CLI and output to stdout
    """
    cmd = ["gbtask", "modify"] + args + ["-i", infile]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
    return result.stdout


@pytest.mark.parametrize(
    ("args", "infile", "expect_outfile"),
    [
        # Test several options of modify.
        (
            [
                "--pad-left",
                "10",
                "--pad-right",
                "xyz",
                "-a",
                "source",
                "test",
                "--annotation",
                "organism",
                "mus musculus",
                "-q",
                "exon",
                "gene",
                "genequalifierforexons",
                "-q",
                "all",
                "locus_tag",
                "locustagqualifierforall",
            ],
            *case_io("test01"),
        ),
        # test --include
        (["--include", "exon"], *case_io("test02")),
        # test --exclude
        (["--exclude", "exon"], *case_io("test03")),
        # test --feature
        (
            [
                "--feature",
                "0",
                "10",
                "1",
                "exon",
                "--feature",
                "-0",
                "-10",
                "-1",
                "exon",
            ],
            *case_io("test04"),
        ),
    ],
)
def test_modify_integration(args, infile, expect_outfile):
    result = modify_run(args, infile)
    with open(expect_outfile, "r") as ef:
        assert result == ef.read()


if __name__ == "__main__":
    # Generate output for test cases
    # Write results of running integration tests to expected outfiles
    integration.generate_integration_output(
        test_modify_integration,
        modify_run,
        infile="infile",
        expect_outfile="expect_outfile",
        process_infile=lambda inf: integration.roundtrip(inf, "genbank"),
        process_outfile=lambda outf: integration.roundtrip(outf, "genbank"),
    )
