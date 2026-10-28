"""Command-line entry point for region annotation."""

import argparse

from regionreport.analysis import annotate_regions
from regionreport.io import read_bed
from regionreport.report import write_report


def main(argv=None):
    """Run the CLI; invalid input is reported by argparse with exit status 2."""
    parser = argparse.ArgumentParser(prog="regionreport", description=__doc__)
    parser.add_argument("command", choices=["annotate"])
    parser.add_argument("--regions", required=True, help="Regions BED4 file")
    parser.add_argument("--genes", required=True, help="Genes BED4 file")
    parser.add_argument("--repeats", required=True, help="TE annotations BED4 file")
    parser.add_argument("--out", required=True, help="New or empty output directory")
    parser.add_argument("--no-plots", action="store_true", help="Write CSV files without PNGs")
    args = parser.parse_args(argv)
    try:
        regions = read_bed(args.regions)
        genes = read_bed(args.genes)
        repeats = read_bed(args.repeats)
        result = annotate_regions(regions, genes, repeats)
        make_plots = not args.no_plots
        directory = write_report(result, args.out, plots=make_plots)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Annotated {len(result.summary)} regions. Report: {directory.resolve()}")
    return 0
