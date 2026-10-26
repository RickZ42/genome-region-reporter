Genome Region Reporter
======================

A Python API and CLI for describing overlaps between genomic regions, genes
and transposable elements. It does not infer biological effects or mechanisms.

Quick start
-----------

Install the package from the repository root::

   python -m pip install ".[docs]"

Run the bundled synthetic example::

   regionreport annotate --regions examples/regions.bed --genes examples/genes.bed --repeats examples/repeats.bed --out results/

Inputs must be tab-separated BED4 files without a header. Coordinates are
zero-based and half-open. Names must be unique within each input; genes have
one interval per gene. Empty inputs are valid. Chromosomes are matched exactly,
and all inputs must describe the same genome assembly. Overlaps ignore strand.

Outputs are ``summary.csv``, ``overlaps.csv``, ``plots.csv`` and optional PNGs.
The output directory must be new or empty. Plot files have numeric names,
mapped to the original region IDs in ``plots.csv``. Summing pairwise overlap
lengths does not give unique coverage when features overlap each other.

API reference
-------------

.. automodule:: regionreport.io
   :members: read_bed, validate_bed

.. automodule:: regionreport.analysis
   :members: AnnotationResult, annotate_regions

.. automodule:: regionreport.report
   :members: write_report

.. automodule:: regionreport.plotting
   :members: plot_region
