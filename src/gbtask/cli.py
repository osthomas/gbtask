import argparse

from gbtask import __version__
from gbtask.geneify import cli as geneify
from gbtask.modify import cli as modify
from gbtask.splice import cli as splice
from gbtask.tofasta import cli as tofasta
from gbtask.togxf import cli as togxf


class Formatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        """
        Limit help width in wide terminals for readability.
        """
        import shutil

        width = shutil.get_terminal_size().columns
        width = min(100, width)
        super().__init__(width=width, *args, **kwargs)

    def _get_help_string(self, action):
        """
        The argparse.ArgumentDefaultsHelpFormatter checks for the presence of
        the formatting string '%(default)' to detect whether a default help
        should be added.
        This variant checks for 'default: ' to allow manual specification of
        a default.
        """
        from gettext import gettext as _

        help = action.help
        if help is None:
            help = ""

        if "default: " not in help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += _(" (default: %(default)s)")
        return help


def parse_args():
    parser = argparse.ArgumentParser(prog="gbtask", formatter_class=Formatter)

    parser.add_argument(
        "-v",
        action="count",
        default=0,
        help="Increase log level. Can be repeated.",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
        help="Print version and exit",
    )
    import sys

    io_parser = argparse.ArgumentParser(description="Input/Output", add_help=False)
    io_parser.add_argument(
        "-i",
        "--infile",
        help="Input file (default: stdin)",
        type=argparse.FileType("r"),
        default=sys.stdin,
    )
    io_parser.add_argument(
        "-o",
        "--outfile",
        help="Output file (default: stdout)",
        type=argparse.FileType("w"),
        default=sys.stdout,
    )
    subparsers = parser.add_subparsers(title="Commands", required=True)

    geneify.get_subparser(subparsers, formatter_class=Formatter, parents=[io_parser])
    modify.get_subparser(subparsers, formatter_class=Formatter, parents=[io_parser])
    splice.get_subparser(subparsers, formatter_class=Formatter, parents=[io_parser])

    togxf.get_subparser_togff(
        subparsers, formatter_class=Formatter, parents=[io_parser]
    )
    togxf.get_subparser_togtf(
        subparsers, formatter_class=Formatter, parents=[io_parser]
    )

    tofasta.get_subparser_tofasta(
        subparsers, formatter_class=Formatter, parents=[io_parser]
    )

    return parser.parse_args()
