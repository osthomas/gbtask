import subprocess
from importlib.resources import files

import pytest

from tests import integration


def case_io(prefix: str):
    """
    Get input and output files for a test case, assuming a common prefix.
    """
    import tests

    inf = files(tests).joinpath(
        "resources", "integration_geneify", prefix + "_input.gb"
    )
    outf = files(tests).joinpath(
        "resources", "integration_geneify", prefix + "_output.gb"
    )
    return inf, outf


def geneify_run(args, infile):
    """
    Run the geneify CLI and output to stdout
    """
    cmd = ["gbtask", "geneify"] + args + ["-i", infile]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
    return result.stdout


@pytest.mark.parametrize(
    ("args", "infile", "expect_outfile"),
    [
        # Group by overlap
        (["--groupby-overlap"], *case_io("test01")),
        # Default Parameters
        ([], *case_io("test02")),
        # Group by nonexistent qualifier, should result in generation of new
        # root features
        (["--groupby", "xyz"], *case_io("test03")),
        # Infer exons and UTRs
        (["--infer-exons", "--infer-utrs"], *case_io("test04")),
    ],
)
def test_geneify_integration(args, infile, expect_outfile):
    result = geneify_run(args, infile)
    with open(expect_outfile, "r") as ef:
        assert result == ef.read()


if __name__ == "__main__":
    # Write results of running integration tests to expected outfiles
    integration.generate_integration_output(
        test_geneify_integration,
        geneify_run,
        infile="infile",
        expect_outfile="expect_outfile",
        process_infile=lambda inf: integration.roundtrip(inf, "genbank"),
        process_outfile=lambda outf: integration.roundtrip(outf, "genbank"),
    )
