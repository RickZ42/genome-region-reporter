"""Pure interval annotation; no files are read or written here."""

from dataclasses import dataclass

import bioframe as bf
import pandas as pd

from regionreport.io import validate_bed

OVERLAP_COLUMNS = [
    "region_id",
    "chrom",
    "region_start",
    "region_end",
    "feature_type",
    "feature_id",
    "feature_start",
    "feature_end",
    "overlap_start",
    "overlap_end",
    "overlap_bp",
]


@dataclass
class AnnotationResult:
    """Summary rows in input region order, plus one row per overlapping pair."""

    summary: pd.DataFrame
    overlaps: pd.DataFrame


def _find_overlaps(regions, features, feature_type):
    if regions.empty or features.empty:
        return pd.DataFrame(columns=OVERLAP_COLUMNS)
    pairs = bf.overlap(
        regions,
        features,
        how="inner",
        return_overlap=True,
        suffixes=("_region", "_feature"),
    )
    return pd.DataFrame(
        {
            "region_id": pairs["name_region"],
            "chrom": pairs["chrom_region"],
            "region_start": pairs["start_region"],
            "region_end": pairs["end_region"],
            "feature_type": feature_type,
            "feature_id": pairs["name_feature"],
            "feature_start": pairs["start_feature"],
            "feature_end": pairs["end_feature"],
            "overlap_start": pairs["overlap_start"],
            "overlap_end": pairs["overlap_end"],
            "overlap_bp": pairs["overlap_end"] - pairs["overlap_start"],
        }
    )


def annotate_regions(
    regions: pd.DataFrame,
    genes: pd.DataFrame,
    repeats: pd.DataFrame,
) -> AnnotationResult:
    """Find positive-length, strand-independent overlaps on matching chromosomes.

    Each argument is a BED4 DataFrame. Counts refer to distinct feature IDs
    within each annotation file. Adjacent half-open intervals do not overlap.
    Names are matched exactly: ``chr1`` and ``1`` are different chromosomes.
    All input tables must use the same genome assembly. Inputs are not mutated.
    """
    regions, genes, repeats = (validate_bed(table) for table in (regions, genes, repeats))
    tables = [
        _find_overlaps(regions, genes, "gene"),
        _find_overlaps(regions, repeats, "TE"),
    ]
    nonempty = [table for table in tables if not table.empty]
    overlaps = (
        pd.concat(nonempty, ignore_index=True)
        if nonempty
        else pd.DataFrame(columns=OVERLAP_COLUMNS)
    )
    if not overlaps.empty:
        order = {name: index for index, name in enumerate(regions["name"])}
        overlaps = (
            overlaps.assign(_order=overlaps["region_id"].map(order))
            .sort_values(["_order", "feature_type", "feature_start", "feature_id"], kind="stable")
            .drop(columns="_order")
            .reset_index(drop=True)
        )

    summary = regions.rename(columns={"name": "region_id"}).copy()
    summary["length_bp"] = summary["end"] - summary["start"]
    for feature_type, column in [("gene", "gene_count"), ("TE", "te_count")]:
        counts = (
            overlaps.loc[overlaps["feature_type"] == feature_type]
            .groupby("region_id")["feature_id"]
            .nunique()
        )
        summary[column] = summary["region_id"].map(counts).fillna(0).astype("int64")
    return AnnotationResult(summary=summary, overlaps=overlaps)
