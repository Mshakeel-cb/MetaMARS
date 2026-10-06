# metaMARS preparation (`metamars_prep`)

A direct Python pipeline for preparing metagenomes and isolate genomes for downstream AMR evidence layers. The Python package has no runtime Python dependencies. Bioinformatics executables and reference databases are installed separately.

This first module produces contigs, genome bins, quality reports, taxonomy, and comparison groups. AMR annotation, transfer inference, and population statistics belong to subsequent evidence modules.

## Workflow

```text
Metagenomic paired reads → fastp → metaSPAdes → canonical contigs
    → read mapping / coverage → MetaBAT2 → CheckM2
    → eligible genomes → GTDB-Tk → taxonomic groups
    + MetaQUAST assembly evaluation

Isolate paired reads → fastp → SPAdes --isolate → canonical genome
    → QUAST + CheckM2 → eligible genomes → GTDB-Tk → groups

Supplied isolate genomes or MAGs → canonical genomes → QUAST + CheckM2
    → eligible genomes → GTDB-Tk → groups
```

Every assembled contig is retained, including short and unbinned sequences and sequences from genomes excluded by quality thresholds. A sample yielding no bins still produces its contig evidence tables. See [the evidence model](docs/evidence-model.md).

## Install and run

Requires Python 3.10 or newer. Run from this directory without installing the Python package:

```bash
python3 -m metamars_prep --help
```

The executable `bin/metamars_prep` also works directly from a source checkout. To install the `metamars_prep` console command in your active Python environment:

```bash
python3 -m pip install .
metamars_prep --help
```

The CLI runs preparation directly. A sample sheet is required for one or multiple inputs.

```text
  -h, --help            help message
  -o, --outdir          output directory
  -t, --threads         working threads dflt=2
  -j, --jobs            Maximum concurrent samples dflt=1
  -m, --memory          memory limit dflt=8
  -tr, --tax-rank       automatically group all classified genomes at specified rank {domain, phylum, class,
                        order,family,genus,species} dflt=species
  -st, --skip-taxonomy  Explicitly omit GTDB-Tk, groups are marked unresolved dflt=False
  --gtdbtk-db           GTDB-Tk reference database directory
  -it, --input-type      {metagenome_reads,isolate_reads,isolate_genomes,MAGs}
  -s, --samples         CSV/TSV sample sheet; paths are relative to this file
```

Accepted input types are exactly `metagenome_reads`, `isolate_reads`, `isolate_genomes`, and `MAGs`. Taxonomic ranks are `domain`, `phylum`, `class`, `order`, `family`, `genus`, and `species`.

## External tools and databases

Example Conda specifications are provided in `envs/`. They use version constraints, not solved lockfiles; validate environment resolution and biological runs on your target machine.

```bash
conda env create -f envs/preparation.yml
conda env create -f envs/checkm2.yml
conda env create -f envs/gtdbtk.yml
```

metaMARS first uses executables available on `PATH`. For missing commands, it automatically looks for the named `metamars-preparation`, `metamars-checkm2`, and `metamars-gtdbtk` environments and launches the tool through `conda run` or `mamba run`. Conda or Mamba must be available on `PATH` (or through `CONDA_EXE`/`MAMBA_EXE`). Keep the names in the supplied YAML files. Environments and databases are never installed or downloaded during a pipeline run.

The environments separate the preparation tools from the CheckM2 and GTDB-Tk dependency stacks. SPAdes includes metaSPAdes, QUAST includes MetaQUAST, and MetaBAT2 includes `jgi_summarize_bam_contig_depths`. Installing the Python package alone does not install these programs.

CheckM2 needs a compatible DIAMOND database, supplied through `CHECKM2DB`. GTDB-Tk uses `--gtdbtk-db` or `GTDBTK_DATA_PATH`; the command-line path takes precedence. GTDB-Tk and its database are unnecessary when `--skip-taxonomy` is set. Databases are not downloaded automatically. Follow the [CheckM2 database instructions](https://github.com/chklovski/CheckM2#database) and [GTDB-Tk installation instructions](https://ecogenomics.github.io/GTDBTk/installing/index.html).

```bash
export CHECKM2DB=/databases/checkm2/uniref100.KO.1.dmnd
export GTDBTK_DATA_PATH=/databases/gtdb
```

## Sample sheets

Put read and genome paths in the sample sheet. Paths are resolved relative to that file. FASTQ and FASTA files may be plain or gzip-compressed.

For mixed input types, provide `input_type` on every row and omit `-it`:

```csv
sample_id,input_type,read1,read2,genome,subject_id,cohort
meta01,metagenome_reads,reads/meta01_R1.fastq.gz,reads/meta01_R2.fastq.gz,,person01,case
iso01,isolate_reads,reads/iso01_R1.fastq.gz,reads/iso01_R2.fastq.gz,,person02,control
iso02,isolate_genomes,,,genomes/iso02.fna,person03,control
mag01,MAGs,,,genomes/mag01.fna,person01,case
```

```bash
metamars_prep -s samples.csv -o results -t 16 -j 2 -m 64 -tr genus
```

For a homogeneous sheet, `-it` supplies the type when the `input_type` column is absent or a cell is blank. Every nonblank `input_type` cell must agree with `-it`.

```csv
sample_id,read1,read2,cohort
sample01,reads/sample01_R1.fastq.gz,reads/sample01_R2.fastq.gz,case
sample02,reads/sample02_R1.fastq.gz,reads/sample02_R2.fastq.gz,control
```

```bash
metamars_prep -s reads.csv -it metagenome_reads -o results -t 8 -m 32
```

One input uses the same format with one data row:

```csv
sample_id,genome
isolate01,genomes/isolate01.fna
```

```bash
metamars_prep -s isolate.csv -it isolate_genomes -o isolate-results \
  --gtdbtk-db /databases/gtdb
```

Read inputs require separate `read1` and `read2` files. Supplied `isolate_genomes` and `MAGs` use the `genome` field and skip raw-read QC, assembly, mapping, and binning. Each row represents one independent input unit; each supplied genome file represents one genome.

Sample IDs must be unique, at most 128 characters, and contain only letters, numbers, underscores, periods or hyphens, starting with a letter or number. Additional columns are preserved as `metadata.*` columns in `samples.tsv`. Use distinct input IDs and a shared `biological_sample_id` for multiple isolates from the same biological sample.

Use an output path containing only letters, numbers, `/`, `.`, `_` and `-`, because external tools impose path restrictions. The mixed sheet in [examples/samples.template.csv](examples/samples.template.csv) is a template to adapt. [examples/samples.csv](examples/samples.csv) illustrates one supplied genome using a tiny synthetic FASTA; it is not a biological benchmark.

## Processing behavior

Memory is specified as an integer number of GB (default 8). The total thread and memory budgets are shared across concurrent sample workers. SPAdes receives the memory budget; for other tools it is advisory, not an operating-system limit. Tools can create auxiliary I/O threads, and GTDB-Tk may need substantially more memory than the default budget.

The preparation settings are fixed: fastp qualified-base Phred threshold 20, minimum read length 50, at least one surviving pair, MetaBAT2 minimum contig length 2,500 bp, and comparison eligibility of at least 50% CheckM2 completeness and at most 10% contamination. The fastp threshold defines a qualified base, not a minimum average read quality. Its default unqualified-base percentage rule remains in effect; adapter trimming is enabled and quality-tail clipping is not enabled.

Host removal is not implemented. Automatic reference downloads are disabled. Assembly evaluation uses reference-free metrics. All original contigs and excluded genomes are retained. These comparison gates do not constitute a full MIMAG quality classification.

`--tax-rank genus` automatically groups eligible genomes by classified genus. Genomes unresolved at the selected rank receive no group ID. `--skip-taxonomy` records `taxonomy_skipped` for otherwise eligible genomes. Tool output streams to the terminal with sample/stage prefixes and is saved in logs.

Samples are assembled independently. Coassembly, cross-sample differential coverage, single-end-only or long-read assembly, and automatic lane merging are not implemented. fastp orphan reads are retained for later analysis but are not assembled or mapped. SPAdes graphs and `contigs.paths` are preserved; `original_id` links canonical contigs to original SPAdes names.

## Running and failures

Use a new or empty `--outdir` for each run. metaMARS executes the pipeline from the beginning; there is no resume option, checkpoint store, output hashing, or automatic retry. Existing nonempty output directories are rejected to protect their contents.

Commands and tool messages stream to the terminal and are saved under `logs/`. A failed tool stops the run and cancels other active sample processes. Partial output and logs remain for diagnosis; rerun into a new directory after correcting the error. `summary.json` and `report.html` are produced after the analyses complete.

## Code structure

`cli.py` parses options, `models.py` reads sample sheets, and `pipeline.py` runs the scientific steps. `runner.py` launches subprocesses and captures logs. `tools.py` finds executables and database paths. `registry.py` preserves contig identity and bin membership; `taxonomy.py` reads quality/taxonomy reports and groups genomes. `config.py` holds run settings and fixed QC thresholds.

## Outputs and validation

Start with `report.html`, `summary.json`, and `run_info.json`. `run_info.json` records inputs, settings, command locations and database paths; it is not execution state. See [the output contract](docs/outputs.md) for table definitions and [the evidence model](docs/evidence-model.md) for interpretation. Relative paths in tables refer to the run output directory.

```bash
python3 -m unittest discover -s tests -v
```

The standard-library tests cover input validation, sequence identity, bin membership, retained unbinned contigs, report parsing, grouping, process failures, environment discovery, and CLI flows with synthetic tool fixtures. The recorded validation also includes real preparation tools on small generated data; see [validation.md](docs/validation.md) for results and the remaining full-pipeline validation limits. No full biological benchmark is bundled.
