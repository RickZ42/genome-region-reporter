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

Outputs are ``summary.csv``, ``overlaps.csv``, ``plots.csv``, ``run.json`` and optional PNGs.
The output directory must be new or empty. Plot files have numeric names,
mapped to the original region IDs in ``plots.csv``. Summing pairwise overlap
lengths does not give unique coverage when features overlap each other.
Summary union-coverage columns merge overlapping intervals before measuring bases.

Circular and GFF3 annotations
----------------------------------------

``read_gff3`` reads plain or gzipped annotations, converts coordinates to BED,
and retains gene-level IDs, labels, source feature types and strand. By default
it selects ``gene`` and ``pseudogene`` records, avoiding child CDS/exon double
counting. Multipart selected IDs are explicitly rejected. Gene-span coverage
includes introns, not just coding bases.

``annotate_regions(..., circular=True, sequence_lengths={"seqid": length})``
accepts intervals crossing the coordinate origin, either with ``start > end`` or
an unwrapped ``end > length``. Lengths are mandatory and intervals may span no
more than one lap. Equal endpoints are rejected; ``[0, length)`` means the whole
chromosome. Returned coordinates are unwrapped. Segments from the same feature
count once per region. All sequences within a call use the chosen topology.

``species`` adds a descriptive label only. It never changes coordinate rules.
The optional ``repeats=None`` records that no mobile-element track was supplied.
The result does not infer biological effects, gene disruption or mobility.

The CLI exposes ``--input-format gff3``, ``--circular``, ``--species``,
``--seqid``, ``--sequence-lengths`` and ``--plot-format png|svg|pdf``. Queries
remain BED4; GFF3 sequence-region directives supply lengths, or a two-column TSV
can be supplied. Strand is displayed but does not filter overlaps.
Length files have no header and use positive integers; blank lines and ``#``
comments are allowed. Invalid rows include a filename and line number in the
error message. Conflicting lengths for the same sequence are rejected.

The repository's ``DEMO/README.md`` documents two real-data examples, including
source URLs, checksums, code, full-chromosome summaries and PNG/SVG images.
``plot_overview`` requires contiguous non-overlapping windows on one sequence;
``plot_region`` caps displayed features and labels any omissions.

API reference
-------------

.. automodule:: regionreport.io
   :members: read_bed, validate_bed, read_sequence_lengths

.. automodule:: regionreport.analysis
   :members: AnnotationResult, annotate_regions

.. automodule:: regionreport.gff
   :members: read_gff3

.. automodule:: regionreport.report
   :members: write_report

.. automodule:: regionreport.plotting
   :members: plot_region, plot_overview
