import json

import pandas as pd
import pytest

from regionreport.cli import main
from regionreport.io import read_sequence_lengths


def test_read_lengths_with_comments_and_repeated_matching_id(tmp_path):
    path = tmp_path / "lengths.tsv"
    path.write_text(
        "\ufeff# Sequence lengths\n\nchr1\t100\nchr2 200\nchr1\t100\n", encoding="utf-8"
    )
    assert read_sequence_lengths(path) == {"chr1": 100, "chr2": 200}


@pytest.mark.parametrize(
    "bad_row",
    [
        "chr2",
        "chr2 100 extra",
        "chr2 zero",
        "chr2 1.5",
        "chr2 -1",
        "chr2 +100",
        "chr2 1_000",
        "chr2 １００",
        "chr2 0",
        f"chr2 {2**63}",
        "chr1 200",
    ],
)
def test_bad_length_reports_filename_and_physical_line(tmp_path, bad_row):
    path = tmp_path / "lengths.tsv"
    path.write_text(f"# comment\nchr1\t100\n\n{bad_row}\n", encoding="utf-8")
    with pytest.raises(ValueError) as error:
        read_sequence_lengths(path)
    assert f"{path}:4:" in str(error.value)


def test_empty_length_file_is_rejected(tmp_path):
    path = tmp_path / "empty.tsv"
    path.write_text("# No records\n\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no sequence lengths found"):
        read_sequence_lengths(path)


def test_largest_supported_length(tmp_path):
    path = tmp_path / "lengths.tsv"
    path.write_text(f"chr1\t{2**63 - 1}\n", encoding="utf-8")
    assert read_sequence_lengths(path) == {"chr1": 2**63 - 1}


def test_cli_bad_length_creates_no_report(tmp_path, capsys):
    regions = tmp_path / "regions.bed"
    genes = tmp_path / "genes.bed"
    regions.write_text("chr1\t0\t10\tquery\n", encoding="utf-8")
    genes.write_text("", encoding="utf-8")
    lengths = tmp_path / "lengths.tsv"
    lengths.write_text("chr1\t0\n", encoding="utf-8")
    output = tmp_path / "report"
    with pytest.raises(SystemExit) as error:
        main(
            [
                "annotate",
                "--regions",
                str(regions),
                "--genes",
                str(genes),
                "--sequence-lengths",
                str(lengths),
                "--out",
                str(output),
            ]
        )
    assert error.value.code == 2
    assert f"{lengths}:1:" in capsys.readouterr().err
    assert not output.exists()


def test_cli_circular_query_uses_length_file(tmp_path):
    regions = tmp_path / "regions.bed"
    genes = tmp_path / "genes.bed"
    lengths = tmp_path / "lengths.tsv"
    regions.write_text("c\t90\t10\tquery\n", encoding="utf-8")
    genes.write_text("c\t95\t5\tgene\n", encoding="utf-8")
    lengths.write_text("c\t100\n", encoding="utf-8")
    output = tmp_path / "report"
    assert (
        main(
            [
                "annotate",
                "--regions",
                str(regions),
                "--genes",
                str(genes),
                "--sequence-lengths",
                str(lengths),
                "--circular",
                "--no-plots",
                "--out",
                str(output),
            ]
        )
        == 0
    )
    summary = pd.read_csv(output / "summary.csv")
    assert summary.loc[0, "length_bp"] == 20
    assert summary.loc[0, "gene_count"] == 1
    assert summary.loc[0, "gene_covered_bp"] == 10
    metadata = json.loads((output / "run.json").read_text(encoding="utf-8"))
    assert metadata["sequence_lengths"] == {"c": 100}
