"""Small fixtures with hand-calculated expected overlaps."""

import pandas as pd
import pytest


def make_bed(rows):
    columns = ["chrom", "start", "end", "name"]
    return pd.DataFrame(rows, columns=columns)


@pytest.fixture
def bed():
    """Give tests a function they can call as bed(rows)."""
    return make_bed
