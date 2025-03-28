import argparse

from Bio import SeqIO


def get_subparser_tofasta(subparsers, *args, **kwargs) -> None:
    parser = subparsers.add_parser(
        "tofasta",
        help="Convert sequences in Genbank file to FASTA",
        *args,
        **kwargs,
    )

    parser.set_defaults(func=tofasta)


def tofasta(args: argparse.Namespace):
    for record in SeqIO.parse(args.infile, "genbank"):
        args.outfile.write(record.format("fasta"))
