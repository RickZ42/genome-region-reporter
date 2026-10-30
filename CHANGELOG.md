# Changelog

## 1.6.0 — 30 October 2026

- Sequence-length files now report the filename and physical line number for
  malformed records. Lengths must be positive integers within the supported
  coordinate range; conflicting records and empty length files are rejected.
- Added `--max-features` to control the number of features shown in each plot.
  The default remains 60. CSV output always retains every overlapping feature.
- Invalid plot limits are rejected before any report files are created.
- `run.json` now records the plot limit and whether plots were enabled.
- Added `regionreport --version` and `python -m regionreport --version`.
- Aligned the package and Python API version numbers with this release; previous
  GitHub releases still reported package version `0.1.0` when installed.

The circular-coordinate, GFF3, union-coverage and real-data demo functionality
introduced in the preceding release is retained.

## Earlier releases

- [v1.5](https://github.com/RickZ42/genome-region-reporter/releases/tag/Genome_region_reporter_v1.5):
  circular sequences, GFF3 annotation input and real-data demonstrations.
- [v1.0](https://github.com/RickZ42/genome-region-reporter/releases/tag/v1.0):
  initial BED region reports, tests and documentation.
