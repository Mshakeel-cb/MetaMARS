# Validation of the simplified preparation pipeline

Validated on 2026-09-25 with Python 3.13.5:

```text
PYTHONDONTWRITEBYTECODE=1 python3 -W error::ResourceWarning -m unittest discover -s tests -q
Ran 100 tests in 9.075s
OK
```

The tests cover CLI and sample-sheet validation, paired FASTQ consistency,
compressed inputs, canonical sequence identity, retained unbinned contigs, exact
bin membership, CheckM2 and GTDB report parsing, taxonomic grouping, subprocess
failures/cancellation, Conda/Mamba command discovery, and synthetic CLI runs.

Integration cases exercise metagenomic reads, isolate reads, supplied genomes,
MAGs, mixed inputs, parallel sample processing, one-thread mapping, zero bins,
all-short contigs, quality exclusions, skipped taxonomy, documented taxonomy
failures, and missing classification reports. Existing output directories are
preserved; failed runs retain logs and must be rerun into a fresh directory.
The removed resume option is rejected. There are no checkpoint records.

Regression checks include R1/R2 hardlinks to the same file, malformed quoted
sample sheets and reports, missing coverage values, and Mamba's different
streaming command syntax. Optional absent trailing metadata remains accepted.

## Real preparation-tool smoke test

Existing local tool installations were also executed on generated data: 12,000
paired 150-base reads sampled from twelve 25,000-base contigs. This used the same
QC, mapping, sorting, indexing, depth, binning, and reference-free assembly
assessment arguments as the pipeline. It did not assemble the generated reads.

| Program | Result |
|---|---|
| fastp 1.1.0 | Produced valid cleaned pairs, orphan files and JSON/HTML reports |
| minimap2 2.31-r1302 + samtools | Produced a sorted BAM and separate BAM index |
| jgi_summarize_bam_contig_depths | Produced positive coverage for all 12 contigs |
| MetaBAT2 2.18 | Recovered one 275,000-base bin; membership validation retained all 12 contigs, including the unbinned contig |
| MetaQUAST/QUAST 5.3.0 | Produced reference-free TSV/HTML reports through a real `conda run` invocation |

This smoke test validates actual command execution and output parsing on small
synthetic data, not biological accuracy or performance on metagenomes.

## Remaining validation limits

SPAdes and GTDB-Tk were unavailable in the inspected installations. CheckM2 was
not run with a real reference database. Their command and output contracts are
covered by synthetic fixtures; the command syntax was also checked against the
[SPAdes documentation](https://ablab.github.io/spades/running.html),
[CheckM2 documentation](https://github.com/chklovski/CheckM2), and
[GTDB-Tk 2.7 documentation](https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html).
The assembly-report options were checked against the
[QUAST manual](https://quast.sourceforge.net/docs/manual.html).

The supplied Conda environment specifications were not installed or solved as a
complete set. No full run with all real tools and databases or biological
benchmark was completed. Wheel installation was not tested because the active
Python environment lacks pip/setuptools. The source CLI and executable launcher
work without installing the Python package.
