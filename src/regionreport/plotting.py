"""Static region plots with a separate labelled row for each feature."""

from pathlib import Path

from matplotlib.figure import Figure


def plot_region(result, region_id, path):
    """Save a PNG showing a region and clipped overlapping annotations.

    Plotting uses a Figure directly and requires no graphical desktop. Labels
    longer than 45 characters are abbreviated; CSV files retain full IDs.
    """
    matches = result.summary.loc[result.summary["region_id"] == region_id]
    if matches.empty:
        raise ValueError(f"Unknown region ID: {region_id}")
    region = matches.iloc[0]
    hits = result.overlaps.loc[result.overlaps["region_id"] == region_id]
    row_count = len(hits) + 1
    height = max(3, 1.8 + 0.42 * row_count)
    figure = Figure(figsize=(10, height), layout="constrained")
    axis = figure.subplots()
    start = int(region["start"])
    end = int(region["end"])

    # Row zero shows the query; each following row shows one annotation.
    axis.barh(0, end - start, left=start, height=0.55, color="#375a7f")
    labels = [f"Region: {region_id}"]
    for index, hit in enumerate(hits.itertuples(index=False), start=1):
        colour = "#c16e29"
        if hit.feature_type == "gene":
            colour = "#25866a"
        width = hit.overlap_end - hit.overlap_start
        axis.barh(index, width, left=hit.overlap_start, height=0.55, color=colour)
        labels.append(f"{hit.feature_type}: {hit.feature_id}")

    for index, label in enumerate(labels):
        if len(label) > 45:
            labels[index] = label[:42] + "..."
    axis.set_yticks(range(row_count), labels=labels, parse_math=False)
    axis.invert_yaxis()
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
