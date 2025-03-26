import argparse
import sys

from Bio import SeqIO

from gbtask.splice import splice


def get_subparser(subparsers, *args, **kwargs) -> None:
    parser = subparsers.add_parser(
        "splice",
        help="Remove parts of a Genbank file not covered by an exon. No attempts are made to disentangle genes. If a location is covered by an exon feature, it is retained.",
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
        "-e", "--exon", help="Name of feature type to consider as exon", default="exon"
    )

    parser.set_defaults(func=run)


def run(args: argparse.Namespace):
    for record in SeqIO.parse(args.infile, "genbank"):
        spliced = splice.splice(record, args.exon)
        SeqIO.write(spliced, args.outfile, "genbank")
