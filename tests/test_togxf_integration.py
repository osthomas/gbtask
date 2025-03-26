import subprocess
from importlib.resources import files

import pytest

import tests
from tests import integration


def case_io_gff(prefix: str):
    """
    Get input and output files for a test case, assuming a common prefix.
    """
    inf = files(tests).joinpath("resources", "integration_togxf", prefix + "_input.gb")
    outf = files(tests).joinpath(
        "resources", "integration_togxf", prefix + "_output.gff"
    )
    return inf, outf


def case_io_gtf(prefix: str):
    """
    Get input and output files for a test case, assuming a common prefix.
    """
    inf = files(tests).joinpath("resources", "integration_togxf", prefix + "_input.gb")
    outf = files(tests).joinpath(
        "resources", "integration_togxf", prefix + "_output.gtf"
    )
    return inf, outf


def togff_run(args, infile):
    """
    Run the modify CLI and output to stdout
    """
    cmd = ["gbtask", "togff"] + args + ["-i", infile]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
    return result.stdout


def togtf_run(args, infile):
    """
    Run the modify CLI and output to stdout
    """
    cmd = ["gbtask", "togtf"] + args + ["-i", infile]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
    return result.stdout


@pytest.mark.parametrize(
    ("args", "infile", "expect_outfile"),
    [
        # Test several options of modify.
        (
            [],
            *case_io_gff("test01"),
        ),
    ],
)
def test_togff_integration(args, infile, expect_outfile):
    result = togff_run(args, infile)
    with open(expect_outfile, "r") as ef:
        assert result == ef.read()


@pytest.mark.parametrize(
    ("args", "infile", "expect_outfile"),
    [
        # Test several options of modify.
        (
            [],
            *case_io_gtf("test01"),
        ),
    ],
)
def test_togtf_integration(args, infile, expect_outfile):
    result = togtf_run(args, infile)
    with open(expect_outfile, "r") as ef:
        assert result == ef.read()


if __name__ == "__main__":
    # Generate output for test cases
    # Write results of running integration tests to expected outfiles
    integration.generate_integration_output(
        test_togff_integration,
        togff_run,
        infile="infile",
        expect_outfile="expect_outfile",
        process_infile=lambda inf: integration.roundtrip(inf, "genbank"),
    )
    integration.generate_integration_output(
        test_togtf_integration,
        togtf_run,
        infile="infile",
        expect_outfile="expect_outfile",
        process_infile=lambda inf: integration.roundtrip(inf, "genbank"),
    )
