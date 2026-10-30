"""Pure interval annotation with explicit linear or circular sequence topology."""

import bioframe as bf
import pandas as pd

from regionreport.io import BED_COLUMNS, validate_bed
from regionreport.topology import check_lengths, split_at_origin, unwrap_intervals

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
    "feature_strand",
    "feature_label",
    "annotation_type",
]
INTEGER_COLUMNS = [
    "region_start",
    "region_end",
    "feature_start",
    "feature_end",
    "overlap_start",
    "overlap_end",
    "overlap_bp",
]


class AnnotationResult:
    """Keep summary, per-segment overlaps and run metadata together.

    Circular overlap rows use an unwrapped coordinate axis relative to the query.
    One feature can contribute multiple segments, but is counted only once per
    region and feature type. No inference of mobility or biological effect occurs.
    """

    def __init__(self, summary, overlaps, metadata=None):
        self.summary = summary
        self.overlaps = overlaps
        self.metadata = {} if metadata is None else metadata


def _empty_overlaps():
    result = pd.DataFrame(columns=OVERLAP_COLUMNS)
    for column in INTEGER_COLUMNS:
        result[column] = result[column].astype("int64")
    return result


def _prepare_features(features, regions, circular, lengths):
    if not circular:
        pieces = features.copy()
        pieces["original_start"] = pieces["start"]
        pieces["original_end"] = pieces["end"]
        return pieces
    # Splitting is essential: duplicating whole wrapped features can lose their
    # low-coordinate portion when the query itself does not cross the origin.
    pieces = split_at_origin(features, lengths)
    shifted = []
    for chrom, queries in regions.groupby("chrom", sort=False):
        length = lengths[chrom]
        extent = int(queries["end"].max()) - length
        if extent <= 0:
            continue
        low = pieces[(pieces["chrom"] == chrom) & (pieces["start"] < extent)].copy()
        # Clip before adding the sequence length to avoid int64 overflow for a
        # long feature beyond the queried part of the second lap.
        low["end"] = low["end"].clip(upper=extent)
        low["start"] = low["start"] + length
        low["end"] = low["end"] + length
        shifted.append(low)
    if shifted:
        pieces = pd.concat([pieces, *shifted], ignore_index=True)
    return pieces


def _find_overlaps(regions, features, feature_type, circular, lengths):
    if regions.empty or features.empty:
        return _empty_overlaps()
    features = _prepare_features(features, regions, circular, lengths)
    features = features.copy()
    for column, default in [("strand", "."), ("annotation_type", feature_type)]:
        if column not in features:
            features[column] = default
    if "label" not in features:
        features["label"] = features["name"]
    pairs = bf.overlap(
        regions,
        features,
        how="inner",
        return_overlap=True,
        suffixes=("_region", "_feature"),
    )
    if pairs.empty:
        return _empty_overlaps()
    return pd.DataFrame(
        {
            "region_id": pairs["name_region"],
            "chrom": pairs["chrom_region"],
            "region_start": pairs["start_region"],
            "region_end": pairs["end_region"],
            "feature_type": feature_type,
            "feature_id": pairs["name_feature"],
            "feature_start": pairs["original_start_feature"],
            "feature_end": pairs["original_end_feature"],
            "overlap_start": pairs["overlap_start"],
            "overlap_end": pairs["overlap_end"],
            "overlap_bp": pairs["overlap_end"] - pairs["overlap_start"],
            "feature_strand": pairs["strand_feature"],
            "feature_label": pairs["label_feature"],
            "annotation_type": pairs["annotation_type_feature"],
        }
    )


def _union_length(starts, ends):
    total, current_end = 0, -1
    for start, end in sorted(zip(starts, ends, strict=True)):
        start, end = int(start), int(end)
        total += max(0, end - max(start, current_end))
        current_end = max(current_end, end)
    return total


def annotate_regions(
    regions,
    genes,
    repeats=None,
    *,
    circular=False,
    sequence_lengths=None,
    species=None,
):
    """Annotate BED tables using explicit topology and unique-feature counts.

    ``circular=True`` requires a sequence-length mapping for every input sequence.
    Circular [start, end) can wrap (start > end) or use an unwrapped end > length.
    Intervals may span at most one lap. Equal endpoints are rejected: use [0, L)
    for a full chromosome. All sequences in a call share the chosen topology;
    process linear and circular replicons separately when topology is mixed.

    ``species`` is descriptive metadata only, never a topology switch. ``repeats``
    may be None when no mobile-element track was supplied; metadata distinguishes
    that absence from an explicitly supplied empty track. Gene/TE covered bases
    use interval unions. Inputs are copied and never modified.
    """
    if not isinstance(circular, bool):
        raise ValueError("circular must be True or False")
    if species is not None and (not isinstance(species, str) or not species.strip()):
        raise ValueError("species must be a nonempty string or None")
    lengths = check_lengths(sequence_lengths)
    if circular and not lengths:
        raise ValueError("circular=True requires sequence_lengths")
    repeats_supplied = repeats is not None
    if repeats is None:
        repeats = pd.DataFrame(columns=BED_COLUMNS)
    frames = []
    for frame in (regions, genes, repeats):
        checked = validate_bed(frame, circular=circular)
        frames.append(unwrap_intervals(checked, circular, lengths))
    regions, genes, repeats = frames
    tables = [
        _find_overlaps(regions, repeats, "TE", circular, lengths),
        _find_overlaps(regions, genes, "gene", circular, lengths),
    ]
    nonempty = [table for table in tables if not table.empty]
    overlaps = pd.concat(nonempty, ignore_index=True) if nonempty else _empty_overlaps()
    summary = regions.rename(columns={"name": "region_id"})
    summary["length_bp"] = summary["end"] - summary["start"]
    for kind, prefix in [("gene", "gene"), ("TE", "te")]:
        hits = overlaps[overlaps["feature_type"] == kind]
        counts = hits.groupby("region_id", sort=False)["feature_id"].nunique()
        coverage = {
            region_id: _union_length(group["overlap_start"], group["overlap_end"])
            for region_id, group in hits.groupby("region_id", sort=False)
        }
        summary[f"{prefix}_count"] = summary["region_id"].map(counts).fillna(0).astype("int64")
        covered = summary["region_id"].map(coverage).fillna(0).astype("int64")
        summary[f"{prefix}_covered_bp"] = covered
        summary[f"{prefix}_coverage_fraction"] = covered / summary["length_bp"]
    if not overlaps.empty:
        order = {name: i for i, name in enumerate(regions["name"])}
        overlaps = (
            overlaps.assign(
                _order=overlaps["region_id"].map(order),
                _type=overlaps["feature_type"].map({"TE": 0, "gene": 1}),
            )
            .sort_values(["_order", "_type", "feature_start", "feature_id", "overlap_start"])
            .drop(columns=["_order", "_type"])
            .reset_index(drop=True)
        )
    metadata = {
        "species": species,
        "circular": circular,
        "sequence_lengths": lengths,
        "coordinate_system": "0-based half-open; circular ends unwrapped past sequence length",
        "mobile_annotations_supplied": repeats_supplied,
        "gene_records": len(genes),
        "mobile_records": len(repeats),
    }
    return AnnotationResult(summary, overlaps, metadata)
