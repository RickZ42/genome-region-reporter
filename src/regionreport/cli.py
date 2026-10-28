"""Command-line entry point for region annotation."""

import argparse
from pathlib import Path

from regionreport.analysis import annotate_regions
from regionreport.io import read_bed
from regionreport.report import write_report


def main(argv: list[str] | None = None) -> int:
    """Run the CLI; invalid input is reported by argparse with exit status 2."""
    parser = argparse.ArgumentParser(prog="regionreport", description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    annotate = subparsers.add_parser("annotate", help="Annotate regions from three BED4 files")
    for name in ["regions", "genes", "repeats"]:
        annotate.add_argument(f"--{name}", type=Path, required=True, help=f"{name} BED4 file")
    annotate.add_argument("--out", type=Path, required=True, help="New or empty output directory")
    annotate.add_argument("--no-plots", action="store_true", help="Write CSV files without PNGs")
    args = parser.parse_args(argv)
    try:
        result = annotate_regions(
            read_bed(args.regions),
            read_bed(args.genes),
            read_bed(args.repeats),
        )
        directory = write_report(result, args.out, plots=not args.no_plots)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Annotated {len(result.summary)} regions. Report: {directory.resolve()}")
    return 0
