import argparse

from Bio import SeqIO

from gbtask.togxf.GXFConverter import GFFConverter, GTFConverter, GXFConverter


def get_subparser_togff(subparsers, *args, **kwargs) -> None:
    parser = subparsers.add_parser(
        "togff",
        help="Convert genbank to gff. Requires proper qualifiers to be set, run through 'geneify' first!",
        *args,
        **kwargs,
    )

    parser.set_defaults(func=togff)


def get_subparser_togtf(subparsers, *args, **kwargs) -> None:
    parser = subparsers.add_parser(
        "togtf",
        help="Convert genbank to gtf. Requires proper qualifiers to be set, run through 'geneify' first!",
        *args,
        **kwargs,
    )

    parser.set_defaults(func=togtf)


def togxf(args: argparse.Namespace, gxfconv=GXFConverter):
    for record in SeqIO.parse(args.infile, "genbank"):
        converter = gxfconv(record)
        for line in converter.convert():
            args.outfile.write(line)
            args.outfile.write("\n")


def togff(args: argparse.Namespace):
    togxf(args, GFFConverter)


def togtf(args: argparse.Namespace):
    togxf(args, GTFConverter)
