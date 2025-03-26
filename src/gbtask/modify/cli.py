import argparse
import sys

from Bio import SeqIO
from gbtask.modify.modify import modify


def get_subparser(subparsers, *args, **kwargs) -> None:
    parser = subparsers.add_parser(
        "modify",
        help="Perform simple modifications on a genbank file",
        *args,
        **kwargs,
    )
    parser.add_argument(
        "-i",
        "--infile",
        help="Input file",
        type=argparse.FileType("r"),
        default=sys.stdin,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Output file",
        type=argparse.FileType("w"),
        default=sys.stdout,
    )

    parser.add_argument(
        "-a",
        "--annotation",
        action="append",
        metavar=("key", "value"),
        nargs=2,
        default=[],
        help="Set annotation 'key' to 'value' for all records. "
        "Can be given multiple times. "
        "Example: -a organism 'mus musculus'",
    )

    parser.add_argument(
        "-q",
        "--qualifier",
        action="append",
        metavar=("type", "key", "value"),
        nargs=3,
        default=[],
        help="For all features of type 'type', set qualifier 'key' to 'value'. "
        "Can be given multiple times. "
        "The special value 'all' for 'type' can be given to affect all features."
        "Example: -q all locus_tag geneA",
    )

    parser.add_argument(
        "-f",
        "--feature",
        action="append",
        metavar=("start", "end", "strand", "type"),
        nargs=4,
        default=[],
        help="Add a new feature. Can be given multiple times. "
        "NOTE: Locations 'start' and 'end' are 0-based and end-exclusive! "
        "If they start with '-', counting begins from the end of the sequence. "
        "'0' refers to the first, '-0' refers to the last base of the record.",
    )

    grp_inex_ = parser.add_argument_group(
        "Include/Exclude Features",
        description="Feature types to include or exclude. "
        "Only --include or --exclude can be specified, not both. "
        "By default, all feature types are included.",
    )
    grp_inex = grp_inex_.add_mutually_exclusive_group()
    grp_inex.add_argument(
        "--include", nargs="*", default=[], help="Only output these feature types."
    )
    grp_inex.add_argument(
        "--exclude",
        nargs="*",
        default=[],
        help="Exclude these feature types from the output.",
    )
    grp_pad = parser.add_argument_group("Padding")
    grp_pad.add_argument(
        "--pad-left",
        help="Pad sequence to the left. This can either be a number to add this many 'N' characters, or a string to pad with a specific sequence.",
    )
    grp_pad.add_argument(
        "--pad-right",
        help="Like --pad-left, but for the right side of the sequence.",
    )
    parser.set_defaults(func=run)


def run(args: argparse.Namespace):
    for record in SeqIO.parse(args.infile, "genbank"):
        modified = modify(
            record,
            pad_left=args.pad_left,
            pad_right=args.pad_right,
            include=args.include,
            exclude=args.exclude,
            features=args.feature,
            annotations=args.annotation,
            qualifiers=args.qualifier,
        )
        SeqIO.write(modified, args.outfile, "genbank")
