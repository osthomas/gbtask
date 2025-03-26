import logging
from importlib.resources import Package, files
from typing import Callable

from Bio import SeqIO

log = logging.getLogger(__name__)


def case_files(
    package: Package,
    case: str,
    infile_suffix: str = "_input.gb",
    outfile_suffix: str = "_output.gb",
):
    """
    Return paths to input file and expected output for a test case.
    """
    dir = files(package)
    infile = str(dir.joinpath(case + infile_suffix))
    expect_file = str(dir.joinpath(case + outfile_suffix))
    return infile, expect_file


def roundtrip(fpath: str, format: str):
    """
    Roundtrip through Bio.SeqIO to avoid test failures from minor formatting
    differences after manual preparation of input files..
    """
    sr = SeqIO.read(fpath, format)
    sr.annotations["date"] = "01-JAN-2000"
    sr.name = "TESTCASE"
    try:
        del sr.annotations["references"]
    except KeyError:
        pass
    sr.features = [f for f in sr.features if f.type != "source"]
    SeqIO.write(sr, fpath, format)


def generate_integration_output(
    func_marked: Callable,
    func_generate: Callable,
    infile: str = "infile",
    expect_outfile: str = "expect_outfile",
    process_infile: Callable = lambda inf: True,
    process_outfile: Callable = lambda outf: True,
):
    """
    Generate output for integration test cases.

    Parameters
    ---------

    func_marked
        The test function marked by @pytest.mark.parametrize
    func_generate
        The function used to generate output for one case of `func_marked`.
        This function receives all arguments specified in @pytest.mark.parametrize
        to `func_marked` **except** the last one.
        The last argument is assumed to be the path to the file containing the
        expected output.
        When called with these arguments, the function must return the output
        expected for `func_marked` when called with the same arguments.
        The return value is written to the file specified by the last argument.
    infile, expect_outfile
        Name of the argument specifying the input and expected output file,
        respectively, as determined in the call to @pytest.mark.parametrize.
        This is used to identify the relevant arguments in the argument list
        for each parametrization.
    process_infile, process_outfile
        Functions that take the path to the input or expected output file,
        respectively, and perform modifications on it, eg. normalization,
        stripping of white space, ...
        The files are expected to be modified in place.
    """
    log.debug(f"Generating output for integration tests from {func_marked.__name__}")
    if "pytestmark" not in func_marked.__dir__():
        raise AttributeError("`func_marked` has no 'pytestmark' attribute")
    mark = None
    for mark in func_marked.pytestmark:  # type: ignore
        if mark.name == "parametrize":
            break
    if mark is None:
        raise AttributeError("`func_marked` is not parametrized")
    arg_names = mark.args[0]
    idx_outf = arg_names.index(expect_outfile)
    idx_inf = arg_names.index(infile)
    for args in mark.args[1]:
        args = list(args)
        outf = args[idx_outf]
        args = args[:idx_outf] + args[idx_outf + 1 :]

        log.debug(f"Processing {args[idx_inf]}")
        process_infile(args[idx_inf])
        log.debug(f"Calling function with {args}")
        result = func_generate(*args)
        with open(outf, "w") as ef:
            ef.write(result)

        log.debug(f"Processing {outf}")
        process_outfile(outf)
