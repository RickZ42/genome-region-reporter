"""Static region plots with a separate labelled row for each feature."""

from pathlib import Path

from matplotlib.figure import Figure

from regionreport.analysis import AnnotationResult


def plot_region(result: AnnotationResult, region_id: str, path: str | Path) -> Path:
    """Save a PNG showing a region and clipped overlapping annotations.

    Plotting uses a Figure directly and requires no graphical desktop. Labels
    longer than 45 characters are abbreviated; CSV files retain full IDs.
    """
    matches = result.summary.loc[result.summary["region_id"] == region_id]
    if matches.empty:
        raise ValueError(f"Unknown region ID: {region_id}")
    region = matches.iloc[0]
    hits = result.overlaps.loc[result.overlaps["region_id"] == region_id]
    rows = [(f"Region: {region_id}", region["start"], region["end"], "#375a7f")]
    for hit in hits.itertuples(index=False):
        colour = "#25866a" if hit.feature_type == "gene" else "#c16e29"
        rows.append(
            (
                f"{hit.feature_type}: {hit.feature_id}",
                hit.overlap_start,
                hit.overlap_end,
                colour,
            )
        )

    figure = Figure(figsize=(10, max(3, 1.8 + 0.42 * len(rows))), layout="constrained")
    axis = figure.subplots()
    for index, (_, start, end, colour) in enumerate(rows):
        axis.barh(index, end - start, left=start, height=0.55, color=colour)
    labels = [label if len(label) <= 45 else label[:42] + "..." for label, *_ in rows]
    axis.set_yticks(range(len(rows)), labels=labels, parse_math=False)
    axis.invert_yaxis()
    start, end = int(region["start"]), int(region["end"])
    padding = max(1, (end - start) * 0.03)
    axis.set_xlim(start - padding, end + padding)
    axis.ticklabel_format(axis="x", style="plain", useOffset=False)
    axis.set_xlabel("Genomic position (bp; 0-based, half-open)")
    axis.set_title(f"{region_id} | {region['chrom']}:[{start}, {end})", parse_math=False)
    axis.grid(axis="x", alpha=0.2)
    path = Path(path)
    try:
        figure.savefig(path, dpi=150, format="png")
    finally:
        figure.clear()
    return path
