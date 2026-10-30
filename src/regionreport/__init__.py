"""Annotate genomic intervals using BED coordinates: zero-based, half-open."""

from regionreport.analysis import AnnotationResult, annotate_regions
from regionreport.gff import read_gff3
from regionreport.io import read_bed
from regionreport.plotting import plot_overview, plot_region
from regionreport.report import write_report

__all__ = [
    "AnnotationResult",
    "annotate_regions",
    "read_bed",
    "read_gff3",
    "write_report",
    "plot_region",
    "plot_overview",
]
__version__ = "0.1.0"
