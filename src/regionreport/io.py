"""Read and validate BED4 tables without changing the caller's data."""

from numbers import Integral
from pathlib import Path

import pandas as pd

BED_COLUMNS = ["chrom", "start", "end", "name"]
MAX_COORDINATE = 2**63 - 1


def read_sequence_lengths(path):
    """Read sequence IDs and positive lengths from a headerless two-column file."""
    path = Path(path)
    lengths = {}
    with path.open(encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            location = f"{path}:{line_number}"
            fields = text.split()
            if len(fields) != 2:
                raise ValueError(f"{location}: expected two columns: sequence ID and length")
            chrom, length_text = fields
            if not length_text.isascii() or not length_text.isdigit():
                raise ValueError(f"{location}: sequence length must be a positive integer")
            length = int(length_text)
            if length < 1 or length > MAX_COORDINATE:
                raise ValueError(
                    f"{location}: sequence length must be between 1 and {MAX_COORDINATE}"
                )
            if chrom in lengths and lengths[chrom] != length:
                raise ValueError(f"{location}: conflicting sequence length for {chrom}")
            lengths[chrom] = length
    if not lengths:
        raise ValueError(f"{path}: no sequence lengths found")
    return lengths


def read_bed(path, *, circular=False):
    """Read the first four BED columns into a checked pandas DataFrame.

    Blank lines, comments and browser settings are skipped. Empty files are valid.
    circular=True permits wrapped start > end; annotate_regions checks lengths.
    """
    path = Path(path)
    rows = []

    with path.open(encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()

            # Skip lines that do not describe a genomic interval.
            if text == "":
                continue
            if text.startswith("#"):
                continue

            first_word = text.split()[0]
            if first_word in ["track", "browser"] and "\t" not in text:
                continue

            # Split on tabs, keeping only chromosome, start, end and name.
            fields = line.rstrip("\r\n").split("\t")
            location = f"{path}:{line_number}"
            if len(fields) < 4:
                raise ValueError(f"{location}: expected at least four tab-separated BED fields")

            chrom = fields[0].strip()
            start_text = fields[1].strip()
            end_text = fields[2].strip()
            name = fields[3].strip()

            # Accept ordinary digits only, then convert them to integers.
            for coordinate in [start_text, end_text]:
                if not coordinate.isascii() or not coordinate.isdigit():
                    raise ValueError(f"{location}: BED coordinates must be nonnegative integers")

            start = int(start_text)
            end = int(end_text)
            invalid_order = start == end if circular else start >= end
            if invalid_order or max(start, end) > MAX_COORDINATE:
                contract = "distinct nonnegative endpoints" if circular else "0 <= start < end"
                raise ValueError(f"{location}: require {contract}, both <= {MAX_COORDINATE}")

            rows.append((chrom, start, end, name))

    frame = pd.DataFrame(rows, columns=BED_COLUMNS)
    try:
        return validate_bed(frame, circular=circular)
    except ValueError as error:
        raise ValueError(f"{path}: {error}") from error


def validate_bed(frame, *, circular=False):
    """Check BED columns and values, returning a copy with integer coordinates."""
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("BED data must be a pandas DataFrame")

    if not frame.columns.is_unique:
        raise ValueError("BED columns must have unique names")
    for column in BED_COLUMNS:
        if column not in frame.columns:
            raise ValueError("BED data needs chrom, start, end and name columns")

    # Work on a copy so the original table is not changed.
    optional = [col for col in ["strand", "label", "annotation_type"] if col in frame.columns]
    result = frame[BED_COLUMNS + optional].copy()
    result = result.reset_index(drop=True)

    # Check text fields and remove spaces at their edges.
    for column in ["chrom", "name"]:
        cleaned_values = []
        for value in result[column]:
            if not isinstance(value, str):
                raise ValueError(f"BED {column} values must be nonempty strings")
            value = value.strip()
            if value == "":
                raise ValueError(f"BED {column} values must be nonempty strings")
            cleaned_values.append(value)
        result[column] = cleaned_values

    if result["name"].duplicated().any():
        raise ValueError("BED names must be unique within each input file")

    # True and False must not be treated as coordinates, even though bool is an int subtype.
    for column in ["start", "end"]:
        for value in result[column]:
            if isinstance(value, bool):
                raise ValueError("BED coordinates must be integers")
            if not isinstance(value, Integral):
                raise ValueError("BED coordinates must be integers")
            if value < 0 or value > MAX_COORDINATE:
                raise ValueError("BED coordinates must be nonnegative signed 64-bit integers")
        result[column] = result[column].astype("int64")

    invalid = result["start"] == result["end"] if circular else result["start"] >= result["end"]
    if invalid.any():
        message = "distinct endpoints in circular mode" if circular else "start < end"
        raise ValueError(f"BED intervals must have {message}")

    for column in optional:
        if not result[column].map(lambda value: isinstance(value, str)).all():
            raise ValueError(f"BED {column} values must be strings")
    if "strand" in result and not result["strand"].isin(["+", "-", ".", "?"]).all():
        raise ValueError("BED strand must be +, -, . or ?")

    # Keep consistent text-column types, including when the table is empty.
    result["chrom"] = result["chrom"].astype("object")
    result["name"] = result["name"].astype("object")
    return result
