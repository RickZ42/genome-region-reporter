import random

import pandas as pd
import pytest

from regionreport import annotate_regions


def test_wrapped_query_and_wrapped_feature_count_once(bed):
    result = annotate_regions(
        bed([("c", 90, 10, "query")]),
        bed([("c", 95, 5, "wrap"), ("c", 0, 8, "low"), ("c", 20, 30, "away")]),
        circular=True,
        sequence_lengths={"c": 100},
        species="Example species",
    )
    row = result.summary.iloc[0]
    assert row["length_bp"] == 20
    assert row["end"] == 110
    assert row["gene_count"] == 2
    assert row["gene_covered_bp"] == 13
    assert row["gene_coverage_fraction"] == 13 / 20
    wrapped = result.overlaps[result.overlaps["feature_id"] == "wrap"]
    assert len(wrapped) == 2
    assert wrapped["overlap_bp"].sum() == 10
    assert wrapped["feature_end"].tolist() == [105, 105]
    assert result.metadata["species"] == "Example species"


def test_nonwrapped_query_sees_low_piece_of_wrapped_feature(bed):
    result = annotate_regions(
        bed([("c", 0, 8, "query")]),
        bed([("c", 95, 105, "gene")]),
        circular=True,
        sequence_lengths={"c": 100},
    )
    assert result.summary.loc[0, "gene_count"] == 1
    assert result.overlaps["overlap_bp"].tolist() == [5]


def test_full_chromosome_does_not_double_count_wrapped_feature(bed):
    result = annotate_regions(
        bed([("c", 0, 100, "whole")]),
        bed([("c", 95, 5, "gene")]),
        circular=True,
        sequence_lengths={"c": 100},
    )
    assert result.summary.loc[0, "gene_count"] == 1
    assert result.summary.loc[0, "gene_covered_bp"] == 10
    assert result.overlaps["overlap_bp"].sum() == 10


def test_union_coverage_not_pairwise_sum(bed):
    result = annotate_regions(
        bed([("c", 0, 100, "q")]),
        bed([("c", 0, 60, "a"), ("c", 40, 100, "b")]),
    )
    assert result.overlaps["overlap_bp"].sum() == 120
    assert result.summary.loc[0, "gene_covered_bp"] == 100
    assert result.summary.loc[0, "gene_coverage_fraction"] == 1


def test_circular_requires_lengths_and_species_is_metadata_only(bed):
    frame = bed([("c", 90, 10, "q")])
    with pytest.raises(ValueError, match="sequence_lengths"):
        annotate_regions(frame, bed([]), circular=True)
    with pytest.raises(ValueError, match="start < end"):
        annotate_regions(frame, bed([]), species="Escherichia coli")
    with pytest.raises(ValueError, match="Missing sequence length"):
        annotate_regions(frame, bed([]), circular=True, sequence_lengths={"other": 100})


@pytest.mark.parametrize("start,end", [(0, 101), (100, 110), (90, 191), (10, 10)])
def test_invalid_circular_intervals(bed, start, end):
    with pytest.raises(ValueError):
        annotate_regions(
            bed([("c", start, end, "q")]), bed([]), circular=True, sequence_lengths={"c": 100}
        )


@pytest.mark.parametrize("length", [0, -1, True, 10.5, 2**63])
def test_invalid_lengths(bed, length):
    with pytest.raises(ValueError):
        annotate_regions(bed([]), bed([]), circular=True, sequence_lengths={"c": length})


def test_linear_length_bound(bed):
    with pytest.raises(ValueError, match="linear sequence length"):
        annotate_regions(bed([("c", 0, 101, "q")]), bed([]), sequence_lengths={"c": 100})


def test_missing_mobile_track_is_distinct_from_empty_track(bed):
    q = bed([("c", 0, 10, "q")])
    assert not annotate_regions(q, bed([])).metadata["mobile_annotations_supplied"]
    assert annotate_regions(q, bed([]), bed([])).metadata["mobile_annotations_supplied"]


def test_random_circular_overlaps_against_basewise_oracle(bed):
    """Independent oracle: explicit sets of bases, without splitting or bioframe."""
    rng = random.Random(314159)
    length = 101

    def make_rows(count, prefix):
        rows = []
        for i in range(count):
            start = rng.randrange(length)
            span = rng.randrange(1, length + 1)
            end = start + span
            if end > length and span < length and rng.choice([True, False]):
                end %= length
            rows.append((rng.choice(["a", "b"]), start, end, f"{prefix}{i}"))
        return rows

    def bases(row):
        _, start, end, _ = row
        if end < start:
            end += length
        return {pos % length for pos in range(start, end)}

    queries = make_rows(35, "q")
    genes, mobile = make_rows(45, "g"), make_rows(25, "m")
    original = bed(queries)
    before = original.copy(deep=True)
    result = annotate_regions(
        original,
        bed(genes),
        bed(mobile),
        circular=True,
        sequence_lengths={"a": length, "b": length},
    )
    pd.testing.assert_frame_equal(original, before)
    for i, query in enumerate(queries):
        q_bases = bases(query)
        assert result.summary.loc[i, "length_bp"] == len(q_bases)
        for kind, prefix, features in [("gene", "gene", genes), ("TE", "te", mobile)]:
            expected = {}
            covered = set()
            for feature in features:
                overlap = q_bases & bases(feature) if query[0] == feature[0] else set()
                if overlap:
                    expected[feature[3]] = len(overlap)
                    covered |= overlap
            actual = (
                result.overlaps[
                    (result.overlaps["region_id"] == query[3])
                    & (result.overlaps["feature_type"] == kind)
                ]
                .groupby("feature_id")["overlap_bp"]
                .sum()
                .to_dict()
            )
            assert actual == expected
            assert result.summary.loc[i, f"{prefix}_count"] == len(expected)
            assert result.summary.loc[i, f"{prefix}_covered_bp"] == len(covered)
