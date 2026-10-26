import subprocess
import sys
from pathlib import Path
from sysconfig import get_path

import pandas as pd
import pytest

from regionreport import annotate_regions, write_report
from regionreport.cli import main
from regionreport.plotting import plot_region


@pytest.fixture
def inputs(tmp_path):
    paths = {}
    for name, text in {
        "regions": "chr1\t1000\t2000\tregion_1\nchr2\t0\t10\tno_hits\n",
        "genes": "chr1\t1200\t1800\tGene_A\n",
        "repeats": "chr1\t1500\t1700\tTE_1\n",
    }.items():
        path = tmp_path / f"{name}.bed"
        path.write_text(text, encoding="utf-8")
        paths[name] = path
    return paths


def arguments(inputs, output):
    return [
        "annotate",
        "--regions",
        str(inputs["regions"]),
        "--genes",
        str(inputs["genes"]),
        "--repeats",
        str(inputs["repeats"]),
        "--out",
        str(output),
    ]


def test_installed_cli_end_to_end(tmp_path, inputs):
    output = tmp_path / "report"
    completed = subprocess.run(
        [sys.executable, "-m", "regionreport", *arguments(inputs, output)],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "Annotated 2 regions" in completed.stdout
    summary = pd.read_csv(output / "summary.csv")
    assert summary["gene_count"].tolist() == [1, 0]
    assert summary["te_count"].tolist() == [1, 0]
    assert sorted(pd.read_csv(output / "overlaps.csv")["overlap_bp"]) == [200, 600]
    plots = pd.read_csv(output / "plots.csv")
    assert plots["region_id"].tolist() == ["region_1", "no_hits"]
    for filename in plots["plot_file"]:
        assert (output / filename).read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_no_plots(tmp_path, inputs):
    output = tmp_path / "csv_only"
    assert main([*arguments(inputs, output), "--no-plots"]) == 0
    assert not (output / "plots").exists()
    assert pd.read_csv(output / "plots.csv").empty


def test_bad_input_creates_no_report(tmp_path, inputs, capsys):
    inputs["genes"].write_text("chr1\t10\t5\tbad\n", encoding="utf-8")
    output = tmp_path / "report"
    with pytest.raises(SystemExit) as error:
        main(arguments(inputs, output))
    assert error.value.code == 2
    assert "start < end" in capsys.readouterr().err
    assert not output.exists()


def test_existing_output_is_protected(tmp_path, inputs):
    output = tmp_path / "report"
    output.mkdir()
    marker = output / "important.txt"
    marker.write_text("keep", encoding="utf-8")
    with pytest.raises(SystemExit) as error:
        main(arguments(inputs, output))
    assert error.value.code == 2
    assert marker.read_text(encoding="utf-8") == "keep"


def test_safe_plot_filenames_and_literal_labels(tmp_path, bed):
    result = annotate_regions(bed([("chr1", 0, 10, "../../$unclosed")]), bed([]), bed([]))
    output = tmp_path / "safe"
    write_report(result, output)
    manifest = pd.read_csv(output / "plots.csv")
    assert manifest.loc[0, "plot_file"] == "plots/region_0001.png"
    assert manifest.loc[0, "region_id"] == "../../$unclosed"
    assert (output / "plots/region_0001.png").is_file()


def test_empty_report(tmp_path, bed):
    output = tmp_path / "empty"
    write_report(annotate_regions(bed([]), bed([]), bed([])), output)
    assert pd.read_csv(output / "summary.csv").empty
    assert pd.read_csv(output / "overlaps.csv").empty
    assert pd.read_csv(output / "plots.csv").empty


def test_plot_unknown_region(tmp_path, bed):
    result = annotate_regions(bed([("chr1", 0, 10, "known")]), bed([]), bed([]))
    with pytest.raises(ValueError, match="Unknown region"):
        plot_region(result, "unknown", tmp_path / "missing.png")


def test_console_entry_point():
    executable = Path(get_path("scripts")) / (
        "regionreport.exe" if sys.platform == "win32" else "regionreport"
    )
    completed = subprocess.run([str(executable), "--help"], capture_output=True, check=True)
    assert b"annotate" in completed.stdout
