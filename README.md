# Genome Region Reporter

A small Python library and command-line tool for annotating genomic regions
with overlapping genes and transposable elements (TEs).

This is an initial teaching implementation for the Oxford Software Engineering
group project. It describes positional overlaps, not biological effects or
the causes of structural variants. All bundled data are synthetic.

## Install

Requires Python 3.10 or newer. From a clone of this repository:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev,docs]"
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1` instead.
For a normal installation without development tools, use `python -m pip install .`.
The package is not published on PyPI.

## Run the example

From the repository root, with the environment activated:

```bash
regionreport annotate --regions examples/regions.bed --genes examples/genes.bed --repeats examples/repeats.bed --out results/
```

`python -m regionreport annotate ...` is an equivalent entry point. Add
`--no-plots` to generate CSV files without PNGs. Use a new or empty output
directory for each run; existing nonempty directories are never overwritten.

Outputs:

| File | Contents |
| --- | --- |
| `summary.csv` | Each region's length and numbers of overlapping genes and TEs |
| `overlaps.csv` | Each overlapping region/feature pair, coordinates and overlap length |
| `plots.csv` | Mapping from original region IDs to numeric PNG filenames |
| `plots/region_0001.png`, ... | One labelled plot per region, including regions with no hits |

Expected counts for the bundled example:

| Region | Genes | TEs |
| --- | ---: | ---: |
| `region_1` | 1 | 2 |
| `region_2` | 1 | 1 |
| `region_no_hits` | 0 | 0 |

`Gene_A` overlaps `region_1` by 600 bp; `TE_1` overlaps it by 200 bp
and `TE_2` by 100 bp. A feature starting exactly where a region ends is
adjacent, not overlapping.

## Input contract

All three inputs are tab-separated BED files with at least four columns:

```text
chrom    start    end    name
```

The line above describes the fields; actual files have **no header**. For example:

```text
chr1	1000	2000	region_1
```

- Coordinates are **0-based, half-open**: `[start, end)`. Length is `end - start`.
- Coordinates must be integers with `0 <= start < end <= 2**63 - 1`.
- Names must be nonempty and unique within each file. Genes should have one
  interval per gene; this version does not resolve transcripts/exons into genes.
- Empty files are valid. Empty annotations yield zero counts; empty regions
  produce header-only CSV files and no PNGs.
- Blank lines, `#` comments, and standard UCSC `track`/`browser` lines are skipped.
- Fields after the fourth are ignored, including strand. Overlaps are strand-independent.
- Chromosome names are matched exactly (`chr1` differs from `1`). All files must
  refer to the same genome assembly; the tool cannot verify assembly identity.
- BED4 gene and TE records are accepted, not raw GFF3, FASTQ, BAM or VCF files.

Feature counts are counts of annotation IDs, not covered bases. In particular,
do not add overlap lengths across overlapping TEs and interpret the sum as
unique TE coverage. Plots clip features to the query region; full feature
coordinates and IDs remain in the CSV. This version targets small regional
reports; plotting every feature can create tall figures for dense regions.

See the [UCSC BED specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
Interval matching uses [Bioframe](https://bioframe.readthedocs.io/en/latest/guide-intervalops.html).

## Python API

```python
from regionreport import annotate_regions, read_bed, write_report

result = annotate_regions(
    read_bed("examples/regions.bed"),
    read_bed("examples/genes.bed"),
    read_bed("examples/repeats.bed"),
)
print(result.summary)
print(result.overlaps)
write_report(result, "results-api", plots=True)
```

`annotate_regions` accepts BED4 pandas DataFrames without reading or writing
files, and leaves its inputs unchanged. `read_bed` and `write_report` handle I/O.

## Development

```bash
python -m ruff check .
python -m ruff format --check .
python -m pytest --cov=regionreport --cov-report=term-missing
python -m sphinx -W --keep-going -b html docs docs/_build/html
python -m build
```

Tests use hand-calculated interval examples, malformed inputs, empty tables,
no-hit regions, safe output filenames, and the installed CLI. The CI workflow
checks Linux, macOS and Windows on Python 3.12, plus Python 3.10 and 3.14 on
Linux. It also builds documentation and installs a built wheel before testing.

The `src/` layout and running tests outside the checkout help catch packaging
problems that an editable local install could hide.

For team contributions, open an issue, work on a branch, and submit a pull
request describing the change and its validation. See [CONTRIBUTING.md](CONTRIBUTING.md).

## Possible next tasks

- Convert a well-defined subset of GFF3 gene annotations into BED4.
- Add nearest-gene distances with an explicit distance convention.
- Compute unique TE coverage using interval unions.
- Add richer examples and test performance on larger annotation files.

These are extensions, not implemented features. Licensed under the MIT licence.
