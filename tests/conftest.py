"""Small fixtures with hand-calculated expected overlaps."""

import pandas as pd
import pytest


@pytest.fixture
def bed():
    """Construct a BED4 DataFrame from rows."""

    def make_bed(rows):
        return pd.DataFrame(rows, columns=["chrom", "start", "end", "name"])

    return make_bed
