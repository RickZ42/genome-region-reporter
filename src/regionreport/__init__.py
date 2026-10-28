"""Annotate genomic intervals using BED coordinates: zero-based, half-open."""

from regionreport.analysis import AnnotationResult, annotate_regions
from regionreport.io import read_bed
from regionreport.report import write_report

__all__ = ["AnnotationResult", "annotate_regions", "read_bed", "write_report"]
__version__ = "0.1.0"
