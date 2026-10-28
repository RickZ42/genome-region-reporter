"""Write CSV results and optional region plots."""

from pathlib import Path

import pandas as pd

from regionreport.plotting import plot_region


def write_report(result, output_dir, plots=True):
    """Write a report to a new or empty directory and return its path.

    Existing nonempty directories are refused to protect previous results.
    ``plots.csv`` maps safe numeric PNG filenames to the original region IDs.
    """
    output_dir = Path(output_dir)
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
            filename = f"plots/region_{index:04d}.png"
            plot_region(result, region_id, output_dir / filename)
            manifest.append((region_id, filename))
    plot_table = pd.DataFrame(manifest, columns=["region_id", "plot_file"])
    plot_table.to_csv(output_dir / "plots.csv", index=False)
    return output_dir
