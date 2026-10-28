import pandas as pd
import pytest

from regionreport import annotate_regions


@pytest.mark.parametrize(
    "chrom, start, end, expected_bp",
    [
        ("chr1", 1200, 1800, 600),  # Feature contained in region.
        ("chr1", 900, 2100, 1000),  # Region contained in feature.
        ("chr1", 900, 1100, 100),
        ("chr1", 1900, 2100, 100),
        ("chr1", 1000, 2000, 1000),
        ("chr1", 1999, 2000, 1),
        ("chr1", 2000, 2100, 0),  # Touching endpoints do not overlap.
        ("chr1", 900, 1000, 0),
        ("chr2", 1200, 1800, 0),
        ("1", 1200, 1800, 0),
    ],
)
def test_overlap_semantics(bed, chrom, start, end, expected_bp):
    regions = bed([("chr1", 1000, 2000, "region")])
    genes = bed([(chrom, start, end, "gene")])
    repeats = bed([])
    result = annotate_regions(regions, genes, repeats)
    expected_count = 0
    if expected_bp > 0:
        expected_count = 1
    assert result.summary.loc[0, "gene_count"] == expected_count
    assert result.summary.loc[0, "te_count"] == 0
    assert result.overlaps["overlap_bp"].sum() == expected_bp


def test_counts_multiple_features_and_preserves_region_order(bed):
    regions = bed([("chr2", 0, 10, "z"), ("chr1", 0, 10, "a"), ("chr3", 0, 5, "none")])
    genes = bed([("chr1", 1, 5, "gene1"), ("chr1", 2, 8, "gene2")])
    repeats = bed([("chr1", 0, 2, "te1"), ("chr2", 0, 4, "te2")])
    result = annotate_regions(regions, genes, repeats)
    assert result.summary["region_id"].tolist() == ["z", "a", "none"]
    assert result.summary["gene_count"].tolist() == [0, 2, 0]
    assert result.summary["te_count"].tolist() == [1, 1, 0]
    assert result.overlaps.loc[0, "region_id"] == "z"


def test_no_input_mutation_and_same_id_in_different_feature_types(bed):
    regions = bed([("chr1", 0, 10, "region")])
    features = bed([("chr1", 2, 4, "same_id")])
    before = features.copy(deep=True)
    result = annotate_regions(regions, features, features)
    pd.testing.assert_frame_equal(features, before)
    assert result.summary.loc[0, "gene_count"] == 1
    assert result.summary.loc[0, "te_count"] == 1


def test_empty_regions(bed):
    result = annotate_regions(bed([]), bed([("chr1", 0, 10, "gene")]), bed([]))
    assert result.summary.empty
    assert result.overlaps.empty
    assert "gene_count" in result.summary


def test_overlapping_queries_each_receive_the_same_feature(bed):
    result = annotate_regions(
        bed([("chr1", 0, 10, "a"), ("chr1", 5, 15, "b")]),
        bed([("chr1", 6, 8, "gene")]),
        bed([]),
    )
    assert result.summary["gene_count"].tolist() == [1, 1]
    assert result.overlaps["overlap_bp"].tolist() == [2, 2]
