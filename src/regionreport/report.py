"""Write CSV results and optional region plots."""

import json
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from regionreport.plotting import plot_region


def write_report(result, output_dir, plots=True, *, plot_format="png", max_features=60):
    """Write a report to a new or empty directory and return its path.

    Existing nonempty directories are refused to protect previous results.
    ``plots.csv`` maps safe numeric plot filenames to the original region IDs.
    ``max_features`` limits plotted features only; CSV tables remain complete.
    """
    output_dir = Path(output_dir)
    if plot_format not in {"png", "svg", "pdf"}:
        raise ValueError("plot_format must be png, svg or pdf")
    if isinstance(max_features, bool) or not isinstance(max_features, int) or max_features < 1:
        raise ValueError("max_features must be a positive integer")
    if output_dir.exists():
        if not output_dir.is_dir() or any(output_dir.iterdir()):
            raise FileExistsError(f"Output must be a new or empty directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    result.summary.to_csv(output_dir / "summary.csv", index=False)
    result.overlaps.to_csv(output_dir / "overlaps.csv", index=False)

    manifest = []
    if plots and not result.summary.empty:
        plot_dir = output_dir / "plots"
        plot_dir.mkdir()
        for index, region_id in enumerate(result.summary["region_id"], start=1):
            filename = f"plots/region_{index:04d}.{plot_format}"
            plot_region(result, region_id, output_dir / filename, max_features=max_features)
            manifest.append((region_id, filename))
    plot_table = pd.DataFrame(manifest, columns=["region_id", "plot_file"])
    plot_table.to_csv(output_dir / "plots.csv", index=False)
    metadata = dict(getattr(result, "metadata", {}))
    metadata["package_version"] = version("genome-region-reporter")
    metadata["dependencies"] = {
        name: version(name) for name in ["bioframe", "pandas", "matplotlib"]
    }
    metadata["plot_format"] = plot_format
    metadata["plots_enabled"] = bool(plots)
    metadata["max_features"] = max_features
    (output_dir / "run.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return output_dir
