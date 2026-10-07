# Preparation output contract (schema 1.1)

Paths in TSV tables are relative to the output directory. Empty values indicate unavailable or unresolved data, never a fabricated measurement. TSV headers are present even when there are zero genomes/bins.

| Output | Meaning |
|---|---|
| `contigs.tsv` | One canonical record per assembly contig, with original ID, sequence SHA256, length, GC, FASTA path and membership status |
| `contig_bin_membership.tsv` | Contig-to-genome relation; empty genome ID means unresolved host |
| `genomes.tsv` | Genome FASTA locations, origin (MAG/isolate), CheckM2 metrics and comparison eligibility |
| `coverage.tsv` | Mean mapped-read depth per canonical contig, or explicit unavailable status |
| `taxonomy.tsv` | Full GTDB taxonomy, domain-to-species fields, method, warnings, status and failure reason |
| `lineage_groups.tsv` | One row per genome, selected rank, taxon, stable group ID and assignment status |
| `groups/<group_id>.tsv` | Members of each assigned taxonomic group |
| `samples.tsv` | Input types, cleaned reads, BAM/BAI, graph and contig-path locations, assembly report directory and preserved metadata |
| `read_qc.tsv` | Raw/clean pair counts, fastp reports and orphan-read paths |
| `samples/<sample>/contigs.fna` | Every validated assembly contig with a canonical identifier |
| `samples/<sample>/unbinned.fna` | Every contig absent from emitted bins; an empty file is valid when all contigs belong to a genome |
| `catalog/genomes/*.fna` | Named genome FASTAs used for CheckM2 and GTDB-Tk |
| `report.html`, `summary.json` | Human-readable and machine-readable run summaries |
| `run_info.json` | Software version, input records, parameters, executable command prefixes and resolved database paths |
| `logs/*.log`, `run.log` | External command/stdout/stderr transcripts and pipeline milestones |
| `work/` | Original tool outputs, including full quality/taxonomy reports, cleaned reads and read alignments |

## Reference databases

Reference databases are shared inputs, not per-sample analysis outputs. When no existing path is supplied through the CLI or its environment-variable fallback, they are reused or downloaded under `--db-dir` (default `~/databases/metamars`):

```text
<db-dir>/
├── checkm2/       # DIAMOND database obtained through the CheckM2 downloader
└── gtdbtk-r232/   # Extracted GTDB-Tk R232 reference package
```

GTDB-Tk data are downloaded directly in Python with checksum verification and safe archive extraction. `--skip-taxonomy` skips their download. `--checkm2-db` and `--gtdbtk-db` select existing database locations instead; their paths take precedence over `CHECKM2DB` and `GTDBTK_DATA_PATH`. Invalid supplied paths stop the run instead of downloading a replacement.

GTDB download progress is visible in the terminal and run logs. Detailed external-tool output, including CheckM2 download output, is always logged and is also displayed with `-v` or `--verbose`. Successfully prepared databases are reused across runs with new output directories; this does not resume sample-processing stages. Database setup may take substantial time and disk space on the first run. `run_info.json` records the resolved paths used by the analysis.

## Identity

An assembly ID is `<sample_id>__asm1`; contigs use `<assembly_id>__contig000001`. MAG identifiers use `<sample_id>__bin000001`; supplied genomes and assembled isolates use `<sample_id>__genome`. These identifiers are stable for unchanged input order/content, not a promise of cross-reassembly biological identity. Always retain the sequence hash, original ID and run provenance when comparing runs.

The contig collection defines sequence identity. Genome FASTAs are views of that collection, not additional independent observations. Later AMR annotation should register one physical locus per canonical contig/coordinate/strand and aggregate it through membership.

## Interpretation

Comparison eligibility is a quality gate, not a statement that the genome is complete or free of errors. Failed quality thresholds never remove sequence evidence. GTDB failures/filters and unexplained missing outputs are distinguished: documented per-genome failures are retained with reasons; unexplained tool omissions fail the run. Explicitly unclassified summaries remain unresolved.

Groups with assigned taxonomy get stable IDs based on the complete lineage through the requested rank. Unresolved, quality-excluded, skipped, and failed classifications have no group ID. Each run creates a fresh `groups/` directory for the chosen rank. No sample-level abundance, prevalence statistics, strain identity, or transfer direction is inferred from these groups.

The coverage column is mapping depth, not normalized taxonomic abundance. Supplied genomes without reads have no inferred depth. GC percentage is `(G + C) / sequence length`, including ambiguous positions in the denominator.

`assembly_graph` refers to the original SPAdes GFA. `contig_paths` points to SPAdes paths for the original contigs. Use `contigs.tsv.original_id` to translate these names; canonicalization does not relabel graph segments.
