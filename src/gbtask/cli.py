import argparse

from gbtask.geneify import cli as geneify
from gbtask.modify import cli as modify
from gbtask.splice import cli as splice
from gbtask.togxf import cli as togxf


def parse_args():
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(prog="gbtask", formatter_class=formatter_class)

    parser.add_argument(
        "-v",
        action="count",
        default=0,
        help="Increase log level. Can be repeated.",
    )
    parser.add_argument("--version", action="count", help="Print version and exit")
    subparsers = parser.add_subparsers(title="Commands", required=True)

    geneify.get_subparser(subparsers, formatter_class=formatter_class)
    modify.get_subparser(subparsers, formatter_class=formatter_class)
    splice.get_subparser(subparsers, formatter_class=formatter_class)

    togxf.get_subparser_togff(subparsers, formatter_class=formatter_class)
    togxf.get_subparser_togtf(subparsers, formatter_class=formatter_class)

    return parser.parse_args()
