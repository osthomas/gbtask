import subprocess
from importlib.resources import files

import pytest

import tests
from tests import integration


def case_io_splice(prefix: str):
    """
    Get input and output files for a test case, assuming a common prefix.
    """
    inf = files(tests).joinpath("resources", "integration_splice", prefix + "_input.gb")
    outf = files(tests).joinpath(
        "resources", "integration_splice", prefix + "_output.gb"
    )
    return inf, outf


def splice_run(args, infile):
    """
    Run the splice CLI and output to stdout
    """
    cmd = ["gbtask", "splice"] + args + ["-i", infile]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
    return result.stdout


@pytest.mark.parametrize(
    ("args", "infile", "expect_outfile"),
    [
        # Test several options of modify.
        (
            [],
            *case_io_splice("test01"),
        ),
    ],
)
def test_splice_integration(args, infile, expect_outfile):
    result = splice_run(args, infile)
    with open(expect_outfile, "r") as ef:
        assert result == ef.read()


if __name__ == "__main__":
    # Generate output for test cases
    # Write results of running integration tests to expected outfiles
    integration.generate_integration_output(
        test_splice_integration,
        splice_run,
        infile="infile",
        expect_outfile="expect_outfile",
        process_infile=lambda inf: integration.roundtrip(inf, "genbank"),
    )
