import argparse
import sys

from Bio import SeqIO


def get_subparser_tofasta(subparsers, *args, **kwargs) -> None:
    parser = subparsers.add_parser(
        "tofasta",
        help="Convert sequences in Genbank file to FASTA",
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

    parser.set_defaults(func=tofasta)


def tofasta(args: argparse.Namespace):
    for record in SeqIO.parse(args.infile, "genbank"):
        args.outfile.write(record.format("fasta"))
