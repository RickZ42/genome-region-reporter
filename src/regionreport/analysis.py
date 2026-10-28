"""Pure interval annotation; no files are read or written here."""

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


class AnnotationResult:
    """Keep the summary and overlap tables together."""

    def __init__(self, summary, overlaps):
        self.summary = summary
        self.overlaps = overlaps


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


def annotate_regions(regions, genes, repeats):
    """Count gene and TE overlaps using 0-based, half-open BED coordinates.

    Inputs must use the same genome assembly. Strand is ignored. The original
    tables are not changed, and regions with no hits receive zero counts.
    """
    regions = validate_bed(regions)
    genes = validate_bed(genes)
    repeats = validate_bed(repeats)

    gene_hits = _find_overlaps(regions, genes, "gene")
    te_hits = _find_overlaps(regions, repeats, "TE")

    summary = regions.rename(columns={"name": "region_id"})
    summary["length_bp"] = summary["end"] - summary["start"]
    summary["gene_count"] = 0
    summary["te_count"] = 0
    ordered_hits = []

    # Process regions in their input order. Every annotation ID is unique
    # within its file, so the number of matching rows is the feature count.
    for index, region_id in enumerate(regions["name"]):
        region_genes = gene_hits[gene_hits["region_id"] == region_id]
        region_tes = te_hits[te_hits["region_id"] == region_id]
        summary.loc[index, "gene_count"] = len(region_genes)
        summary.loc[index, "te_count"] = len(region_tes)

        for hits in [region_tes, region_genes]:
            if not hits.empty:
                hits = hits.sort_values(["feature_start", "feature_id"])
                ordered_hits.append(hits)

    overlaps = pd.DataFrame(columns=OVERLAP_COLUMNS)
    if ordered_hits:
        overlaps = pd.concat(ordered_hits, ignore_index=True)

    return AnnotationResult(summary, overlaps)
