"""Read gene-level GFF3 annotations, including explicit mobile-element features."""

import gzip
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

from regionreport.io import BED_COLUMNS, MAX_COORDINATE, validate_bed

GENE_TYPES = ("gene", "pseudogene")
MOBILE_TYPES = ("mobile_genetic_element", "transposable_element", "insertion_sequence")


def read_gff3(path, *, feature_types=GENE_TYPES, seqids=None):
    """Convert selected GFF3 records to BED4 with strand and display metadata.

    Coordinates change from 1-based inclusive to 0-based half-open. Select genes,
    not their child CDS/exon records, to avoid double counting. IDs are retained
    separately from labels. Multi-row selected IDs are rejected explicitly:
    silently replacing discontiguous parts with their bounding span is incorrect.
    Gzip files and an embedded FASTA section are supported. Sequence lengths and
    circular sequence IDs are exposed in DataFrame.attrs, never inferred by name.
    """
    path = Path(path)
    if isinstance(feature_types, str):
        raise ValueError("feature_types must be a collection, not a single string")
    selected = set(feature_types)
    if not selected:
        raise ValueError("Select at least one GFF3 feature type")
    if isinstance(seqids, str):
        seqids = {seqids}
    wanted = None if seqids is None else set(seqids)
    rows, lengths, circular_ids = [], {}, set()
    seen, observed_sequences = set(), set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, 1):
            text = line.strip()
            location = f"{path}:{line_number}"
            if text == "##FASTA":
                break
            if text.startswith("##sequence-region "):
                fields = text.split()
                if len(fields) != 4:
                    raise ValueError(f"{location}: malformed sequence-region directive")
                chrom = unquote(fields[1])
                try:
                    start, end = int(fields[2]), int(fields[3])
                except ValueError as error:
                    raise ValueError(f"{location}: invalid sequence-region coordinates") from error
                if start != 1 or not 0 < end <= MAX_COORDINATE:
                    raise ValueError(f"{location}: require full sequence-region beginning at 1")
                if wanted is None or chrom in wanted:
                    observed_sequences.add(chrom)
                    if chrom in lengths and lengths[chrom] != end:
                        raise ValueError(f"{location}: conflicting sequence lengths")
                    lengths[chrom] = end
                continue
            if not text or text.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{location}: expected nine tab-separated GFF3 fields")
            chrom = unquote(fields[0])
            if wanted is not None and chrom not in wanted:
                continue
            observed_sequences.add(chrom)
            attributes = {}
            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attributes[key] = unquote(value)
            if fields[2] == "region" and attributes.get("Is_circular", "").lower() == "true":
                circular_ids.add(chrom)
            if fields[2] not in selected:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise ValueError(f"{location}: GFF3 coordinates must be integers") from error
            if not 1 <= start <= end <= MAX_COORDINATE:
                raise ValueError(f"{location}: require 1 <= start <= end <= {MAX_COORDINATE}")
            feature_id = attributes.get("ID") or attributes.get("locus_tag")
            if not feature_id:
                raise ValueError(f"{location}: selected feature needs ID or locus_tag")
            if feature_id in seen:
                raise ValueError(
                    f"{location}: repeated selected ID {feature_id!r}; multipart unsupported"
                )
            seen.add(feature_id)
            label = (
                attributes.get("mobile_element_type")
                or attributes.get("Name")
                or attributes.get("locus_tag")
                or feature_id
            )
            rows.append((chrom, start - 1, end, feature_id, fields[6], label, fields[2]))
    if wanted is not None and wanted - observed_sequences:
        missing = ", ".join(sorted(wanted - observed_sequences))
        raise ValueError(f"Requested sequence IDs not found in GFF3: {missing}")
    frame = pd.DataFrame(rows, columns=[*BED_COLUMNS, "strand", "label", "annotation_type"])
    frame = validate_bed(frame)
    frame.attrs = {"sequence_lengths": lengths, "circular_sequences": sorted(circular_ids)}
    return frame
