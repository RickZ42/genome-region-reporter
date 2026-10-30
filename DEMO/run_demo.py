#!/usr/bin/env python3
"""Reproduce real-annotation demos. Run from any directory after installing the package."""

import argparse
import gzip
import hashlib
import json
import platform
import tempfile
from collections import Counter
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from urllib.request import urlopen

import matplotlib as mpl
import pandas as pd

import regionreport
from regionreport import (
    annotate_regions,
    plot_overview,
    plot_region,
    read_gff3,
    write_report,
)
from regionreport.gff import MOBILE_TYPES

HERE = Path(__file__).resolve().parent


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def chromosome_subset(content, seqid):
    lines = []
    for line in content.decode("utf-8").splitlines(keepends=True):
        if line.startswith("##FASTA"):
            break
        if line.startswith("##sequence-region "):
            if line.split()[1] == seqid:
                lines.append(line)
        elif line.startswith(("#!", "##gff-version", "##species", seqid + "\t")):
            lines.append(line)
    return "".join(lines).encode("utf-8")


def prepare_source(config, download):
    path = HERE / "data" / config["file"]
    if download or not path.exists():
        with urlopen(config["source_url"], timeout=120) as response:
            content = gzip.decompress(response.read())
        if sha256(content) != config["source_content_sha256"]:
            raise ValueError(
                "Upstream annotation changed; refusing an unreviewed replacement snapshot"
            )
        selected = (
            chromosome_subset(content, config["seqid"])
            if config["key"] == "arabidopsis"
            else content
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(gzip.compress(selected, mtime=0))
    if sha256(gzip.decompress(path.read_bytes())) != config["content_sha256"]:
        raise ValueError(f"Annotation checksum mismatch: {path}")
    return path


def windows(seqid, length, step):
    return pd.DataFrame(
        [
            (seqid, start, min(start + step, length), f"window_{i:04d}")
            for i, start in enumerate(range(0, length, step), start=1)
        ],
        columns=["chrom", "start", "end", "name"],
    )


def direct_union_coverage(frame, length, circular):
    """Independent sorted-interval union of original features on one chromosome."""
    segments = []
    for start, end in zip(frame["start"], frame["end"], strict=True):
        start, end = int(start), int(end)
        if circular and end > length:
            segments.extend([(start, length), (0, end - length)])
        else:
            segments.append((start, end))
    covered, right = 0, 0
    for start, end in sorted(segments):
        if end > right:
            covered += end - max(start, right)
            right = end
    return covered


def run_one(config, source, root):
    seqid = config["seqid"]
    genes = read_gff3(source, seqids=[seqid])
    repeats = (
        read_gff3(source, feature_types=MOBILE_TYPES, seqids=[seqid])
        if config["mobile_track_supplied"]
        else None
    )
    lengths = genes.attrs["sequence_lengths"]
    length = lengths[seqid]
    if (seqid in genes.attrs["circular_sequences"]) != config["circular"]:
        raise ValueError("Configured topology disagrees with the source annotation")
    kwargs = dict(circular=config["circular"], sequence_lengths=lengths, species=config["species"])
    processed = HERE / "data" / "processed" / config["key"]
    processed.mkdir(parents=True, exist_ok=True)
    genes.to_csv(processed / "genes.tsv", sep="\t", index=False)
    genes[["chrom", "start", "end", "name"]].to_csv(
        processed / "genes.bed",
        sep="\t",
        header=False,
        index=False,
    )
    if repeats is not None:
        repeats.to_csv(processed / "mobile_elements.tsv", sep="\t", index=False)
    binned = windows(seqid, length, config["window_bp"])
    binned.to_csv(processed / "windows.bed", sep="\t", header=False, index=False)
    overview = annotate_regions(binned, genes, repeats, **kwargs)
    overview.metadata.update(
        assembly=config["assembly"],
        source_url=config["source_url"],
        annotation_sha256=config["content_sha256"],
    )
    output = root / config["key"]
    write_report(overview, output / "windows", plots=False)
    for ext in ("png", "svg"):
        plot_overview(overview, output / f"overview.{ext}")
    if config["circular"]:
        local_rows = [(seqid, length - 6000, 6000, "origin_crossing")]
        first_mobile = repeats.sort_values("start").iloc[0]
        local_rows.append(
            (
                seqid,
                max(0, int(first_mobile["start"]) - 6000),
                min(length, int(first_mobile["end"]) + 6000),
                "mobile_neighborhood",
            )
        )
    else:
        local_rows = [(seqid, 100000, 145000, "linear_gene_neighborhood")]
    local_queries = pd.DataFrame(local_rows, columns=["chrom", "start", "end", "name"])
    local_queries.to_csv(processed / "local_regions.bed", sep="\t", header=False, index=False)
    local = annotate_regions(local_queries, genes, repeats, **kwargs)
    local.metadata.update(overview.metadata)
    write_report(local, output / "local", plots=False)
    for region_id in local_queries["name"]:
        for ext in ("png", "svg"):
            plot_region(local, region_id, output / f"{region_id}.{ext}")
    full = annotate_regions(windows(seqid, length, length), genes, repeats, **kwargs)
    full_summary = full.summary.iloc[0]
    direct_union = direct_union_coverage(genes, length, config["circular"])
    checks = {
        "all_genes_counted_once_on_full_chromosome": int(full_summary["gene_count"]) == len(genes),
        "all_mobile_records_counted_once_on_full_chromosome": int(full_summary["te_count"])
        == (len(repeats) if repeats is not None else 0),
        "gene_union_matches_independent_merge": int(full_summary["gene_covered_bp"])
        == direct_union,
        "window_union_sums_to_full_chromosome_union": int(overview.summary["gene_covered_bp"].sum())
        == direct_union,
        "windows_tile_complete_chromosome": int(overview.summary["length_bp"].sum()) == length,
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    summary = {
        "species": config["species"],
        "assembly": config["assembly"],
        "sequence": seqid,
        "sequence_length_bp": length,
        "circular": config["circular"],
        "gene_records": len(genes),
        "gene_types": dict(Counter(genes["annotation_type"])),
        "mobile_annotations_supplied": repeats is not None,
        "mobile_records": len(repeats) if repeats is not None else None,
        "gene_union_covered_bp": direct_union,
        "window_bp": config["window_bp"],
        "window_count": len(binned),
        "local_regions": local.summary.to_dict("records"),
        "checks": checks,
    }
    (output / "validation.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        f"{config['key']}: {length:,} bp; {len(genes):,} gene-level records; "
        f"{len(repeats) if repeats is not None else 'no supplied'} mobile-element records"
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, default=HERE / "results", help="New or empty output directory"
    )
    parser.add_argument(
        "--download", action="store_true", help="Re-download and verify pinned sources"
    )
    args = parser.parse_args()
    root = args.out.resolve()
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        parser.error(f"Output must be new or empty: {root}")
    configs = json.loads((HERE / "data/sources.json").read_text())
    sources = [prepare_source(config, args.download) for config in configs]
    root.parent.mkdir(parents=True, exist_ok=True)
    # Assemble the full demo in a sibling directory; failed runs leave no partial
    # final report. Existing report directories are always protected.
    with tempfile.TemporaryDirectory(prefix=".regionreport-demo-", dir=root.parent) as temp:
        staging = Path(temp) / "report"
        staging.mkdir()
        with mpl.rc_context(
            {"svg.fonttype": "none", "font.family": "DejaVu Sans", "font.size": 10}
        ):
            summaries = [
                run_one(config, source, staging)
                for config, source in zip(configs, sources, strict=True)
            ]
        manifest = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "python": platform.python_version(),
            "versions": {
                name: version(name)
                for name in ["genome-region-reporter", "bioframe", "pandas", "matplotlib"]
            },
            "sources": configs,
            "implementation_sha256": {
                source.name: sha256(source.read_bytes())
                for source in sorted(Path(regionreport.__file__).parent.glob("*.py"))
            },
            "demo_script_sha256": sha256(Path(__file__).read_bytes()),
            "results": summaries,
            "interpretation": (
                "Descriptive annotation overlap only; no biological effect or mobility inference."
            ),
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        if root.exists():
            root.rmdir()
        staging.rename(root)
    print(f"Report: {root}")


if __name__ == "__main__":
    main()
