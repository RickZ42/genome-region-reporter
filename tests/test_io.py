import pandas as pd
import pytest

from regionreport.io import read_bed, validate_bed


def test_read_bed_with_comments_and_extra_columns(tmp_path):
    path = tmp_path / "example.bed"
    path.write_text(
        '# comment\ntrack name="example"\nbrowser position chr1:1-20\n\n'
        "chr1\t0\t10\tgene_1\t0\t-\n",
        encoding="utf-8",
    )
    frame = read_bed(path)
    assert frame.to_dict("records") == [{"chrom": "chr1", "start": 0, "end": 10, "name": "gene_1"}]


@pytest.mark.parametrize(
    "text",
    [
        "chr1 0 10 id\n",
        "chr1\t0\t10\n",
        "chr1\t-1\t10\tid\n",
        "chr1\t0.5\t10\tid\n",
        "chr1\t10\t10\tid\n",
        "chr1\t20\t10\tid\n",
        "chr1\t0\t9223372036854775808\tid\n",
        "chr1\t+1\t10\tid\n",
        "chr1\t1_0\t20\tid\n",
        "chr1\t١\t10\tid\n",
        "chr1\t²\t10\tid\n",
    ],
)
def test_bad_bed_reports_line_number(tmp_path, text):
    path = tmp_path / "bad.bed"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match=r"bad.bed:1:"):
        read_bed(path)


def test_empty_bed_is_valid(tmp_path):
    path = tmp_path / "empty.bed"
    path.write_text("# no annotations\n", encoding="utf-8")
    frame = read_bed(path)
    assert frame.empty
    assert list(frame.columns) == ["chrom", "start", "end", "name"]
    assert str(frame["start"].dtype) == "int64"


@pytest.mark.parametrize(
    "row",
    [
        ("chr1", 0.0, 10, "id"),
        ("chr1", True, 10, "id"),
        ("chr1", -1, 10, "id"),
        ("chr1", 10, 0, "id"),
        ("", 0, 10, "id"),
        ("chr1", 0, 10, " "),
        (None, 0, 10, "id"),
    ],
)
def test_invalid_dataframe(bed, row):
    with pytest.raises(ValueError):
        validate_bed(bed([row]))


def test_duplicate_ids_are_rejected(bed):
    with pytest.raises(ValueError, match="unique"):
        validate_bed(bed([("chr1", 0, 10, "id"), ("chr2", 5, 15, "id")]))


def test_names_and_chromosomes_are_trimmed(bed):
    checked = validate_bed(bed([(" chr1 ", 0, 10, " gene_A ")]))
    assert checked.loc[0, "chrom"] == "chr1"
    assert checked.loc[0, "name"] == "gene_A"


def test_duplicate_ids_after_trimming_are_rejected(bed):
    with pytest.raises(ValueError, match="unique"):
        validate_bed(bed([("chr1", 0, 10, "id"), ("chr2", 5, 15, " id ")]))


def test_dataframe_is_copied_and_extra_columns_are_ignored(bed):
    frame = bed([("chr1", 0, 10, "id")]).assign(score=0)
    before = frame.copy(deep=True)
    checked = validate_bed(frame)
    checked.loc[0, "start"] = 2
    pd.testing.assert_frame_equal(frame, before)
    assert "score" not in checked


def test_missing_columns():
    with pytest.raises(ValueError, match="columns"):
        validate_bed(pd.DataFrame({"chrom": ["chr1"]}))


def test_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_bed(tmp_path / "missing.bed")
