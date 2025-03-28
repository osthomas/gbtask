import argparse
import logging

from Bio import SeqIO

from gbtask.geneify.geneify import geneify

log = logging.getLogger(__name__)


def get_subparser(subparsers, *args, **kwargs) -> None:
    parser = subparsers.add_parser(
        "geneify",
        help="Organize gene features hierarchically and infer features",
        *args,
        **kwargs,
    )

    parser.add_argument(
        "--groupby",
        default=["gene"],
        nargs="*",
        help="Names of qualifiers for first-level grouping of features",
    )

    parser.add_argument(
        "--groupby-overlap",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Group features by overlap after grouping by qualifiers? Warning, this may lead to undesirable separation of features if there is no scaffolding super-feature that covers all sub-features.",
    )

    parser.add_argument(
        "--infer-exons",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Infer mRNA/exon features per feature group, depending on what is available",
    )

    parser.add_argument(
        "--infer-utrs",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Infer UTR features per feature group based on exons and CDS",
    )
    parser.set_defaults(func=run)


def run(args: argparse.Namespace):
    for record in SeqIO.parse(args.infile, "genbank"):
        genefied = geneify(
            record,
            group_by=args.groupby,
            byoverlap=args.groupby_overlap,
            infer_exons=args.infer_exons,
            infer_utrs=args.infer_utrs,
        )
        SeqIO.write(genefied, args.outfile, "genbank")
