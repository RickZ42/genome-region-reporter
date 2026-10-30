"""Normalize circular intervals without guessing topology from a species name."""

from collections.abc import Mapping
from numbers import Integral

import pandas as pd

from regionreport.io import MAX_COORDINATE


def check_lengths(sequence_lengths):
    """Validate an optional mapping of sequence IDs to positive lengths."""
    if sequence_lengths is None:
        return {}
    if not isinstance(sequence_lengths, Mapping):
        raise ValueError("sequence_lengths must map sequence IDs to positive integers")
    lengths = {}
    for chrom, length in sequence_lengths.items():
        if not isinstance(chrom, str) or not chrom.strip():
            raise ValueError("Sequence IDs must be nonempty strings")
        if isinstance(length, bool) or not isinstance(length, Integral):
            raise ValueError("Sequence lengths must be positive integers")
        if length <= 0 or length > MAX_COORDINATE:
            raise ValueError("Sequence lengths must fit positive signed 64-bit integers")
        lengths[chrom] = int(length)
    return lengths


def unwrap_intervals(frame, circular, lengths):
    """Represent origin-crossing [start, end) as [start, end + length).

    Circular intervals must start inside the sequence and span at most one lap.
    Equal endpoints are ambiguous and rejected; [0, length) means a full lap.
    Native GFF3 intervals with end > length are also accepted in circular mode.
    """
    frame = frame.copy()
    for chrom, rows in frame.groupby("chrom", sort=False):
        length = lengths.get(chrom)
        if circular and length is None:
            raise ValueError(f"Missing sequence length for circular sequence: {chrom}")
        if length is None:
            continue
        for index in rows.index:
            start = int(frame.at[index, "start"])
            end = int(frame.at[index, "end"])
            if not circular:
                if end > length:
                    raise ValueError(f"Interval exceeds linear sequence length: {chrom}")
                continue
            if start >= length:
                raise ValueError(f"Circular start must be smaller than sequence length: {chrom}")
            if end < start:
                end += length
            if not start < end <= start + length:
                raise ValueError(f"Circular interval must span > 0 and <= one lap: {chrom}")
            if end > MAX_COORDINATE:
                raise ValueError("Unwrapped circular coordinate exceeds signed 64-bit range")
            frame.at[index, "end"] = end
    return frame


def split_at_origin(frame, lengths):
    """Split unwrapped features into canonical pieces, keeping original metadata."""
    pieces = []
    for row in frame.to_dict("records"):
        length = lengths[row["chrom"]]
        start, end = int(row["start"]), int(row["end"])
        first = row.copy()
        first["original_start"], first["original_end"] = start, end
        first["end"] = min(end, length)
        pieces.append(first)
        if end > length:
            second = first.copy()
            second["start"], second["end"] = 0, end - length
            pieces.append(second)
    columns = [*frame.columns, "original_start", "original_end"]
    result = pd.DataFrame(pieces, columns=columns)
    for col in ("start", "end", "original_start", "original_end"):
        result[col] = result[col].astype("int64")
    return result
