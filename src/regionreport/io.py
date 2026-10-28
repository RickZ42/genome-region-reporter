"""Read and validate BED4 tables without changing the caller's data."""

import re
from numbers import Integral
from pathlib import Path

import pandas as pd

BED_COLUMNS = ["chrom", "start", "end", "name"]
MAX_COORDINATE = 2**63 - 1


def validate_bed(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a validated BED4 copy with integer coordinates.

    Require nonempty chromosome names, globally unique nonempty feature IDs,
    and integer coordinates satisfying ``0 <= start < end <= 2**63 - 1``.
    Empty tables are valid. Extra columns are discarded; strand is ignored.
    """
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("BED data must be a pandas DataFrame")
    if not frame.columns.is_unique or not set(BED_COLUMNS).issubset(frame.columns):
        raise ValueError("BED data needs unique chrom, start, end and name columns")

    result = frame[BED_COLUMNS].copy().reset_index(drop=True)
    for column in ["chrom", "name"]:
        if (
            not result[column]
            .map(lambda value: isinstance(value, str) and bool(value.strip()))
            .all()
        ):
            raise ValueError(f"BED {column} values must be nonempty strings")
        result[column] = result[column].str.strip()
    if result["name"].duplicated().any():
        raise ValueError("BED names must be unique within each input file")

    for column in ["start", "end"]:
        for value in result[column]:
            if isinstance(value, bool) or not isinstance(value, Integral):
                raise ValueError("BED coordinates must be integers")
            if not 0 <= value <= MAX_COORDINATE:
                raise ValueError("BED coordinates must be nonnegative signed 64-bit integers")
        result[column] = result[column].astype("int64")
    if (result["start"] >= result["end"]).any():
        raise ValueError(
            "BED intervals must have start < end (zero-length intervals are unsupported)"
        )
    return result.astype({"chrom": "object", "name": "object"})


def read_bed(path: str | Path) -> pd.DataFrame:
    """Read a tab-separated BED file with at least four columns and no header.

    Skip blank lines, ``#`` comments and UCSC track/browser directives. Retain
    only the first four fields. Errors include the filename and line number.
    """
    path = Path(path)
    rows = []
    with path.open(encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            if text.split(maxsplit=1)[0] in {"track", "browser"} and "\t" not in text:
                continue
            fields = line.rstrip("\r\n").split("\t")
            location = f"{path}:{line_number}"
            if len(fields) < 4:
                raise ValueError(f"{location}: expected at least four tab-separated BED fields")
            chrom, start, end, name = (field.strip() for field in fields[:4])
            if not re.fullmatch(r"[0-9]+", start) or not re.fullmatch(r"[0-9]+", end):
                raise ValueError(f"{location}: BED coordinates must be nonnegative integers")
            start, end = int(start), int(end)
            if start >= end or end > MAX_COORDINATE:
                raise ValueError(f"{location}: require 0 <= start < end <= {MAX_COORDINATE}")
            rows.append((chrom, start, end, name))
    try:
        return validate_bed(pd.DataFrame(rows, columns=BED_COLUMNS))
    except ValueError as error:
        raise ValueError(f"{path}: {error}") from error
