"""Bounded region plots and chromosome-wide summaries from computed results."""

from pathlib import Path

from matplotlib import rc_context
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, MaxNLocator

COLOURS = {"region": "#375a7f", "gene": "#25866a", "TE": "#c16e29"}


def _save(figure, path):
    path = Path(path)
    suffix = path.suffix.lstrip(".").lower() or "png"
    if suffix not in {"png", "svg", "pdf"}:
        raise ValueError("Plot format must be PNG, SVG or PDF")
    try:
        with rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
            figure.savefig(path, dpi=150, format=suffix)
    finally:
        figure.clear()
    return path


def plot_region(result, region_id, path, *, max_features=60):
    """Plot clipped annotations, preserving strand and circular origin crossings.

    Multiple overlap segments of one feature share a row. At most max_features
    features are displayed; omitted counts are explicit and CSVs remain complete.
    PNG, SVG and PDF are selected by the filename extension.
    """
    if isinstance(max_features, bool) or not isinstance(max_features, int) or max_features < 1:
        raise ValueError("max_features must be a positive integer")
    matches = result.summary.loc[result.summary["region_id"] == region_id]
    if matches.empty:
        raise ValueError(f"Unknown region ID: {region_id}")
    region = matches.iloc[0]
    hits = result.overlaps.loc[result.overlaps["region_id"] == region_id]
    keys = ["feature_type", "feature_id"]
    groups = list(hits.groupby(keys, sort=False))
    groups.sort(
        key=lambda item: (
            0 if item[0][0] == "TE" else 1,
            int(item[1]["overlap_start"].min()),
            item[0][1],
        )
    )
    omitted = max(0, len(groups) - max_features)
    groups = groups[:max_features]
    row_count = len(groups) + 1
    height = max(3.4, 2.0 + 0.30 * row_count)
    figure = Figure(figsize=(12, height), layout="constrained")
    axis = figure.subplots()
    start, end = int(region["start"]), int(region["end"])
    axis.barh(0, end - start, left=start, height=0.55, color=COLOURS["region"])
    labels = [f"Region: {region_id}"]
    for index, ((kind, feature_id), segments) in enumerate(groups, start=1):
        label = segments.iloc[0].get("feature_label", feature_id)
        strand = segments.iloc[0].get("feature_strand", ".")
        strand_label = f" ({strand})" if strand in {"+", "-"} else ""
        labels.append(f"{kind}: {label}{strand_label}")
        for hit in segments.itertuples(index=False):
            left, right = int(hit.overlap_start), int(hit.overlap_end)
            axis.barh(index, right - left, left=left, height=0.55, color=COLOURS[kind])
            if strand in {"+", "-"}:
                tip, tail = (right, left) if strand == "+" else (left, right)
                axis.annotate(
                    "",
                    xy=(tip, index),
                    xytext=(tail, index),
                    arrowprops={"arrowstyle": "->", "color": "#163c3b", "lw": 0.8},
                )
    labels = [label if len(label) <= 45 else label[:42] + "..." for label in labels]
    axis.set_yticks(range(row_count), labels=labels, parse_math=False)
    axis.invert_yaxis()
    padding = max(1, (end - start) * 0.03)
    axis.set_xlim(start - padding, end + padding)
    axis.ticklabel_format(axis="x", style="plain", useOffset=False)
    metadata = getattr(result, "metadata", {})
    length = metadata.get("sequence_lengths", {}).get(region["chrom"])
    wrapped = bool(metadata.get("circular") and length and end > length)
    if wrapped:
        axis.axvline(length, color="#8e526b", ls="--", lw=1.1)
        axis.text(
            length,
            1.01,
            "sequence coordinate 0",
            transform=axis.get_xaxis_transform(),
            ha="center",
            va="bottom",
            color="#8e526b",
            fontsize=9,
        )
        axis.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value) % length:,}"))
    axis.set_xlabel(
        "Sequence position (bp; 0-based, half-open)"
        + (" | ticks wrap at sequence length" if wrapped else "")
    )
    axis.set_title(
        f"{region_id} | {region['chrom']}:[{start:,}, {end:,})",
        parse_math=False,
        pad=26 if wrapped else 12,
    )
    axis.grid(axis="x", alpha=0.2)
    notes = []
    if omitted:
        notes.append(
            f"Showing {max_features} of {len(groups) + omitted} features; CSV is complete."
        )
    if wrapped:
        notes.append("Query and overlap ends use an unwrapped coordinate axis in the CSV.")
    if not metadata.get("mobile_annotations_supplied", True):
        notes.append("Mobile-element annotations were not supplied; absence is not inferred.")
    if notes:
        figure.supxlabel("\n".join(notes), fontsize=9, color="#555555")
    return _save(figure, path)


def plot_overview(result, path):
    """Plot counts and union coverage across ordered, non-overlapping windows.

    The input must cover one sequence with contiguous non-overlapping windows;
    using full-chromosome windows gives a full-chromosome overview. Feature counts
    refer to overlaps per window, so a boundary-spanning feature can occur twice.
    """
    summary = result.summary.sort_values("start")
    if summary.empty or summary["chrom"].nunique() != 1:
        raise ValueError("Overview needs nonempty windows from exactly one sequence")
    if (summary["start"].iloc[1:].to_numpy() != summary["end"].iloc[:-1].to_numpy()).any():
        raise ValueError("Overview windows must be contiguous and non-overlapping")
    metadata = getattr(result, "metadata", {})
    chrom = summary.iloc[0]["chrom"]
    length = metadata.get("sequence_lengths", {}).get(chrom)
    if length is not None and summary["end"].max() > length:
        raise ValueError("Overview windows must use canonical sequence coordinates")
    scale = 1e6 if int(summary["end"].max()) >= 1_000_000 else 1e3
    units = "Mb" if scale == 1e6 else "kb"
    x = (summary["start"] + summary["length_bp"] / 2) / scale
    widths = summary["length_bp"] / scale
    figure = Figure(figsize=(13, 8), layout="constrained")
    axes = figure.subplots(3, 1, sharex=True, gridspec_kw={"height_ratios": [1, 1, 0.8]})
    axes[0].bar(x, summary["gene_count"], width=widths * 0.92, color=COLOURS["gene"])
    axes[0].set_ylabel("Gene overlaps / window")
    coverage = summary["gene_coverage_fraction"] * 100
    axes[1].plot(x, coverage, color=COLOURS["gene"], lw=1.5)
    axes[1].fill_between(x, coverage, alpha=0.15, color=COLOURS["gene"])
    axes[1].set_ylabel("Gene union coverage (%)")
    axes[1].set_ylim(0, 100)
    supplied = metadata.get("mobile_annotations_supplied", True)
    if supplied:
        axes[2].bar(x, summary["te_count"], width=widths * 0.92, color=COLOURS["TE"])
        axes[2].set_ylabel("Mobile-element\noverlaps / window")
    else:
        axes[2].text(
            0.5,
            0.5,
            "Mobile-element track not supplied in this annotation",
            transform=axes[2].transAxes,
            ha="center",
            va="center",
            color="#806044",
        )
        axes[2].set_yticks([])
    for axis in axes:
        axis.grid(axis="y", alpha=0.18)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].yaxis.set_major_locator(MaxNLocator(integer=True))
    if supplied:
        axes[2].yaxis.set_major_locator(MaxNLocator(integer=True))
    axes[2].set_xlim(summary["start"].min() / scale, summary["end"].max() / scale)
    axes[2].set_xlabel(f"{chrom} position ({units}; 0-based)")
    species = metadata.get("species") or "Genome region report"
    topology = "circular" if metadata.get("circular") else "linear"
    window_bp = int(summary["length_bp"].max())
    figure.suptitle(
        f"{species}\n{chrom} | {topology} | "
        f"{metadata.get('gene_records', '?'):,} gene records | "
        f"{len(summary)} windows (up to {window_bp / 1000:g} kb)",
        fontsize=13,
    )
    figure.supxlabel(
        "Counts are per-window overlaps; genes crossing a window boundary appear in both.\n"
        "Coverage is the union of annotated gene intervals, not summed pairwise lengths.",
        fontsize=9,
        color="#555555",
    )
    return _save(figure, path)
