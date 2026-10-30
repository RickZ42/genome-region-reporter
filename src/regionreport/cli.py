"""Command-line entry point for region annotation."""

import argparse

from regionreport import __version__
from regionreport.analysis import annotate_regions
from regionreport.gff import MOBILE_TYPES, read_gff3
from regionreport.io import read_bed, read_sequence_lengths
from regionreport.report import write_report


def main(argv=None):
    """Run the CLI; invalid input is reported by argparse with exit status 2."""
    parser = argparse.ArgumentParser(prog="regionreport", description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("command", choices=["annotate"])
    parser.add_argument("--regions", required=True, help="Regions BED4 file")
    parser.add_argument("--genes", required=True, help="Genes BED4 or GFF3 file")
    parser.add_argument("--repeats", help="Optional mobile annotations, in the chosen input format")
    parser.add_argument("--out", required=True, help="New or empty output directory")
    parser.add_argument("--no-plots", action="store_true", help="Write CSV files without PNGs")
    parser.add_argument(
        "--input-format",
        choices=["bed", "gff3"],
        default="bed",
        help="Annotation format; --regions always uses BED4",
    )
    parser.add_argument("--circular", action="store_true", help="All input sequences are circular")
    parser.add_argument("--species", help="Descriptive species/strain label; does not set topology")
    parser.add_argument("--sequence-lengths", help="Headerless two-column sequence ID/length TSV")
    parser.add_argument("--seqid", action="append", help="Keep this GFF3 sequence ID (repeatable)")
    parser.add_argument("--plot-format", choices=["png", "svg", "pdf"], default="png")
    parser.add_argument(
        "--max-features",
        type=int,
        default=60,
        help="Maximum features shown per plot (default: 60); CSVs always stay complete",
    )
    args = parser.parse_args(argv)
    if args.max_features < 1:
        parser.error("--max-features must be a positive integer")
    try:
        regions = read_bed(args.regions, circular=args.circular)
        lengths = {}
        if args.input_format == "gff3":
            genes = read_gff3(args.genes, seqids=args.seqid)
            repeats = (
                read_gff3(args.repeats, feature_types=MOBILE_TYPES, seqids=args.seqid)
                if args.repeats
                else None
            )
            lengths.update(genes.attrs.get("sequence_lengths", {}))
            if repeats is not None:
                for chrom, length in repeats.attrs.get("sequence_lengths", {}).items():
                    if chrom in lengths and lengths[chrom] != length:
                        raise ValueError(f"Conflicting annotation sequence lengths: {chrom}")
                    lengths[chrom] = length
        else:
            genes = read_bed(args.genes, circular=args.circular)
            repeats = read_bed(args.repeats, circular=args.circular) if args.repeats else None
        if args.sequence_lengths:
            for chrom, length in read_sequence_lengths(args.sequence_lengths).items():
                if chrom in lengths and lengths[chrom] != length:
                    raise ValueError(
                        f"{args.sequence_lengths}: length for {chrom} "
                        "conflicts with GFF3 annotations"
                    )
                lengths[chrom] = length
        result = annotate_regions(
            regions,
            genes,
            repeats,
            circular=args.circular,
            sequence_lengths=lengths,
            species=args.species,
        )
        make_plots = not args.no_plots
        directory = write_report(
            result,
            args.out,
            plots=make_plots,
            plot_format=args.plot_format,
            max_features=args.max_features,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Annotated {len(result.summary)} regions. Report: {directory.resolve()}")
    return 0
