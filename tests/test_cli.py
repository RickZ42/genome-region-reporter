import json
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path
from sysconfig import get_path

import pandas as pd
import pytest
from matplotlib.figure import Figure

from regionreport import __version__, annotate_regions, write_report
from regionreport.cli import main
from regionreport.plotting import plot_region

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


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
    args = ["annotate"]
    for name in ["regions", "genes", "repeats"]:
        args.append("--" + name)
        args.append(str(inputs[name]))
    args.append("--out")
    args.append(str(output))
    return args


def test_installed_cli_end_to_end(tmp_path, inputs):
    output = tmp_path / "report"
    command = [sys.executable, "-m", "regionreport"]
    command.extend(arguments(inputs, output))
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "Annotated 2 regions" in completed.stdout
    summary = pd.read_csv(output / "summary.csv")
    assert summary["gene_count"].tolist() == [1, 0]
    assert summary["te_count"].tolist() == [1, 0]
    overlaps = pd.read_csv(output / "overlaps.csv")
    assert sorted(overlaps["overlap_bp"]) == [200, 600]
    plots = pd.read_csv(output / "plots.csv")
    assert plots["region_id"].tolist() == ["region_1", "no_hits"]
    for filename in plots["plot_file"]:
        assert (output / filename).read_bytes().startswith(PNG_SIGNATURE)


def test_no_plots(tmp_path, inputs):
    output = tmp_path / "csv_only"
    args = arguments(inputs, output)
    args.append("--no-plots")
    assert main(args) == 0
    assert not (output / "plots").exists()
    assert pd.read_csv(output / "plots.csv").empty
    metadata = json.loads((output / "run.json").read_text(encoding="utf-8"))
    assert metadata["plots_enabled"] is False


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


def test_plot_clips_features_and_shortens_long_labels(tmp_path, bed, monkeypatch):
    regions = bed([("chr1", 100, 200, "region")])
    genes = bed([("chr1", 80, 130, "gene_" + "A" * 60)])
    repeats = bed([("chr1", 180, 240, "TE")])
    result = annotate_regions(regions, genes, repeats)

    def check_figure(figure, path, **kwargs):
        axis = figure.axes[0]
        bars = []
        for bar in axis.patches:
            bars.append((bar.get_x(), bar.get_width()))
        # Query first, then TE and gene. Only the intersecting portions are drawn.
        assert bars == [(100, 100), (180, 20), (100, 30)]
        gene_label = axis.get_yticklabels()[2].get_text()
        assert len(gene_label) == 45
        assert gene_label.endswith("...")

    monkeypatch.setattr(Figure, "savefig", check_figure)
    plot_region(result, "region", tmp_path / "region.png")


def test_console_entry_point():
    command_name = "regionreport"
    if sys.platform == "win32":
        command_name = "regionreport.exe"
    executable = Path(get_path("scripts")) / command_name
    completed = subprocess.run([str(executable), "--help"], capture_output=True, check=True)
    assert b"annotate" in completed.stdout
    completed = subprocess.run(
        [str(executable), "--version"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert completed.stdout.strip() == f"regionreport {__version__}"


def test_module_version_needs_no_input_files(tmp_path):
    completed = subprocess.run(
        [sys.executable, "-m", "regionreport", "--version"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )
    assert completed.stdout.strip() == f"regionreport {__version__}"
    assert list(tmp_path.iterdir()) == []


def test_report_version_matches_installed_package(tmp_path, bed):
    assert version("genome-region-reporter") == __version__
    output = tmp_path / "versioned"
    result = annotate_regions(bed([]), bed([]))
    write_report(result, output, plots=False)
    metadata = json.loads((output / "run.json").read_text(encoding="utf-8"))
    assert metadata["package_version"] == __version__


def test_cli_feature_limit_keeps_complete_csvs(tmp_path, inputs):
    output = tmp_path / "limited"
    args = arguments(inputs, output)
    args.extend(["--max-features", "1", "--plot-format", "svg"])
    assert main(args) == 0
    svg = (output / "plots/region_0001.svg").read_text(encoding="utf-8")
    assert "Showing 1 of 2 features; CSV is complete." in svg
    overlaps = pd.read_csv(output / "overlaps.csv")
    assert set(overlaps["feature_id"]) == {"Gene_A", "TE_1"}
    summary = pd.read_csv(output / "summary.csv")
    assert summary.loc[0, "gene_count"] == 1
    assert summary.loc[0, "te_count"] == 1
    metadata = json.loads((output / "run.json").read_text(encoding="utf-8"))
    assert metadata["max_features"] == 1
    assert metadata["plots_enabled"] is True
    assert metadata["plot_format"] == "svg"


@pytest.mark.parametrize("limit", ["0", "-1", "1.5", "many"])
def test_cli_invalid_feature_limit_creates_no_report(tmp_path, inputs, capsys, limit):
    output = tmp_path / "invalid"
    args = arguments(inputs, output)
    args.extend(["--max-features", limit])
    with pytest.raises(SystemExit) as error:
        main(args)
    assert error.value.code == 2
    assert "--max-features" in capsys.readouterr().err
    assert not output.exists()


@pytest.mark.parametrize("limit", [0, -1, True, 1.5, "3", None])
def test_report_invalid_feature_limit_creates_no_files(tmp_path, bed, limit):
    output = tmp_path / "invalid"
    result = annotate_regions(bed([("c", 0, 10, "query")]), bed([]))
    with pytest.raises(ValueError, match="max_features must be a positive integer"):
        write_report(result, output, max_features=limit)
    assert not output.exists()
