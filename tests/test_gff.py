import gzip
import json

import pytest

from regionreport import read_gff3
from regionreport.cli import main
from regionreport.gff import MOBILE_TYPES


@pytest.fixture
def gff(tmp_path):
    path = tmp_path / "annotation.gff3.gz"
    text = (
        "##gff-version 3\n##sequence-region c 1 100\n"
        "c\tsrc\tregion\t1\t100\t.\t+\t.\tID=c;Is_circular=true\n"
        "c\tsrc\tgene\t96\t105\t.\t-\t.\tID=g1;Name=wrap%20gene\n"
        "c\tsrc\tCDS\t96\t105\t.\t-\t0\tID=cds1;Parent=g1\n"
        "c\tsrc\tgene\t20\t30\t.\t+\t.\tID=g2;Name=wrap%20gene\n"
        "c\tsrc\tmobile_genetic_element\t2\t4\t.\t+\t.\t"
        "ID=is1;mobile_element_type=insertion%20sequence:IS1\n"
        "other\tsrc\tgene\t1\t3\t.\t.\t.\tID=other\n"
        "##FASTA\n>c\nACGT\n"
    )
    with gzip.open(path, "wt") as handle:
        handle.write(text)
    return path


def test_gene_level_conversion_and_metadata(gff):
    frame = read_gff3(gff, seqids=["c"])
    assert frame["name"].tolist() == ["g1", "g2"]
    assert frame["start"].tolist() == [95, 19]
    assert frame["end"].tolist() == [105, 30]
    assert frame["label"].tolist() == ["wrap gene", "wrap gene"]
    assert frame["strand"].tolist() == ["-", "+"]
    assert frame.attrs["sequence_lengths"] == {"c": 100}
    assert frame.attrs["circular_sequences"] == ["c"]
    assert "cds1" not in frame["name"].tolist()


def test_mobile_types_are_explicit(gff):
    mobile = read_gff3(gff, feature_types=MOBILE_TYPES, seqids="c")
    assert mobile["name"].tolist() == ["is1"]
    assert mobile.iloc[0]["annotation_type"] == "mobile_genetic_element"
    assert mobile.iloc[0]["label"] == "insertion sequence:IS1"


@pytest.mark.parametrize(
    "line",
    [
        "c\tsrc\tgene\t0\t5\t.\t+\t.\tID=x\n",
        "c\tsrc\tgene\t1\t5\t.\t+\t.\tName=x\n",
        "c\tsrc\tgene\t1\t5\t.\t+\t.\tID=x\nc\tsrc\tgene\t6\t8\t.\t+\t.\tID=x\n",
        "c\tsrc\tgene\t1\t5\t.\twrong\t.\tID=x\n",
    ],
)
def test_bad_selected_features_rejected(tmp_path, line):
    path = tmp_path / "bad.gff3"
    path.write_text(line)
    with pytest.raises(ValueError):
        read_gff3(path)


def test_circular_gff_cli(tmp_path, gff):
    regions = tmp_path / "regions.bed"
    regions.write_text("c\t90\t10\tquery\n")
    output = tmp_path / "report"
    assert (
        main(
            [
                "annotate",
                "--regions",
                str(regions),
                "--genes",
                str(gff),
                "--repeats",
                str(gff),
                "--input-format",
                "gff3",
                "--seqid",
                "c",
                "--circular",
                "--species",
                "Test species",
                "--plot-format",
                "svg",
                "--out",
                str(output),
            ]
        )
        == 0
    )
    assert "<svg" in (output / "plots/region_0001.svg").read_text()
    metadata = json.loads((output / "run.json").read_text())
    assert metadata["circular"] is True
    assert metadata["species"] == "Test species"
    assert metadata["sequence_lengths"] == {"c": 100}


def test_unknown_sequence_filter_is_rejected(gff):
    with pytest.raises(ValueError, match="not found"):
        read_gff3(gff, seqids=["typo"])


def test_length_file_cannot_override_gff_length(tmp_path, gff, capsys):
    regions = tmp_path / "regions.bed"
    regions.write_text("c\t90\t10\tquery\n", encoding="utf-8")
    lengths = tmp_path / "lengths.tsv"
    lengths.write_text("c\t200\n", encoding="utf-8")
    output = tmp_path / "report"
    with pytest.raises(SystemExit) as error:
        main(
            [
                "annotate",
                "--regions",
                str(regions),
                "--genes",
                str(gff),
                "--input-format",
                "gff3",
                "--seqid",
                "c",
                "--circular",
                "--sequence-lengths",
                str(lengths),
                "--out",
                str(output),
            ]
        )
    assert error.value.code == 2
    assert "conflicts with GFF3" in capsys.readouterr().err
    assert not output.exists()
