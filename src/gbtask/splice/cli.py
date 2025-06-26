import argparse

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
        "-e",
        "--exon",
        help="Comma-separated list of feature types to consider as exonic. Example: --exon \"CDS,3'UTR,5'UTR\"",
        default="exon",
    )

    parser.set_defaults(func=run)


def run(args: argparse.Namespace):
    args.exon = args.exon.split(",")
    for record in SeqIO.parse(args.infile, "genbank"):
        spliced = splice.splice(record, args.exon)
        SeqIO.write(spliced, args.outfile, "genbank")
