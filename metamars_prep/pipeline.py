"""Preparation workflow: sequence identity is independent of genome membership."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
import csv
import gzip
import html
import itertools
import json
import logging
import math
import os
from pathlib import Path
import shutil

from . import __version__
from .registry import build_membership, normalize_fasta, write_tsv
from .runner import Runner
from .taxonomy import RANKS, collect_gtdb_failures, group_genomes, parse_checkm2, parse_gtdb
from .tools import discover_tools, required_tools
from .databases import ensure_checkm2, validate_checkm2
from .gtdb_download import ensure_gtdbtk, validate_gtdbtk

LOG = logging.getLogger("metamars")


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, default=str) + "\n")


def worker_environment(threads: int) -> dict[str, str]:
    return {name: str(threads) for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}


def fastq_records(path: Path):
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="ascii") as stream:
        number = 0
        while True:
            name = stream.readline()
            if not name:
                return
            sequence, plus, quality = (stream.readline().rstrip("\r\n") for _ in range(3))
            number += 1
            if not name.startswith("@") or not name[1:].strip() or not plus.startswith("+") or not sequence or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record {number} in {path}")
            if any(base.upper() not in "ACGTRYSWKMBDHVN" for base in sequence):
                raise ValueError(f"Invalid sequence in FASTQ record {number} in {path}")
            if any(ord(char) < 33 or ord(char) > 126 for char in quality):
                raise ValueError(f"Invalid quality in FASTQ record {number} in {path}")
            yield name[1:].split()[0].removesuffix("/1").removesuffix("/2")


def validate_pairs(read1: Path, read2: Path, minimum: int) -> int:
    count = 0
    for left, right in itertools.zip_longest(fastq_records(read1), fastq_records(read2)):
        count += 1
        if left is None or right is None or left != right:
            raise ValueError(f"Mismatched paired FASTQ records at pair {count}: {read1}, {read2}")
    if count < minimum:
        raise ValueError(f"Only {count} read pairs remain; required at least {minimum}: {read1}")
    return count


def coverage_rows(path: Path, contigs: list[dict], sample_id: str) -> list[dict]:
    known = {row["contig_id"] for row in contigs}
    found = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"contigName", "totalAvgDepth"}.issubset(reader.fieldnames or []):
            raise ValueError(f"Invalid MetaBAT depth table: {path}")
        for row in reader:
            name = row["contigName"]
            if name not in known or name in found:
                raise ValueError(f"Unknown or repeated coverage contig: {name}")
            try:
                depth = float(row["totalAvgDepth"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Missing or invalid coverage for {name} in {path}") from exc
            if not math.isfinite(depth) or depth < 0:
                raise ValueError(f"Invalid coverage for {name}")
            found[name] = depth
    return [{"sample_id": sample_id, "contig_id": name, "mean_depth": found.get(name, ""),
             "status": "measured" if name in found else "unavailable"} for name in sorted(known)]


def process_sample(sample, args, runner: Runner, tools: dict, threads: int, memory: int) -> dict:
    """Each worker owns one sample's stage directories and canonical exports."""
    out = args.outdir
    work = out / "work" / "samples" / sample.sample_id
    data = out / "samples" / sample.sample_id
    data.mkdir(parents=True, exist_ok=True)
    env = worker_environment(threads)

    def tool(name):
        return tools[name]

    def run(name, commands, outputs):
        runner.run(f"{sample.sample_id}.{name}", commands, work / name, outputs=outputs, env=env)

    def align(reference, read1, read2):
        bam = work / "mapping" / "reads.bam"
        mapper = [*tool("minimap2"), "-ax", "sr", "-t", str(max(1, threads - 2)),
                  str(reference), str(read1), str(read2)]
        sorter = [*tool("samtools"), "sort", "-@", "0", "-m", f"{max(64, min(1024, memory * 128))}M",
                  "-o", str(bam), "-"]
        run("mapping", [mapper, sorter], [bam])
        return bam

    LOG.info("[%s] Preparing %s", sample.sample_id, sample.input_type)
    clean_reads = []
    qc = {"sample_id": sample.sample_id, "status": "not_applicable", "raw_pairs": "", "clean_pairs": ""}
    graph = contig_paths = ""
    bam = bai = None
    if sample.input_type.endswith("_reads"):
        raw_pairs = validate_pairs(sample.read1, sample.read2, 1)
        directory = work / "qc"
        read1, read2 = directory / "clean.R1.fastq.gz", directory / "clean.R2.fastq.gz"
        unpaired1, unpaired2 = directory / "unpaired.R1.fastq.gz", directory / "unpaired.R2.fastq.gz"
        json_report, html_report = directory / "fastp.json", directory / "fastp.html"
        command = [*tool("fastp"), "--in1", str(sample.read1), "--in2", str(sample.read2),
                   "--out1", str(read1), "--out2", str(read2), "--unpaired1", str(unpaired1), "--unpaired2", str(unpaired2),
                   "--json", str(json_report), "--html", str(html_report), "--thread", str(min(16, threads)),
                   "--detect_adapter_for_pe", "--length_required", str(args.min_read_length),
                   "--qualified_quality_phred", str(args.qualified_phred)]
        run("qc", [command], [read1, read2, unpaired1, unpaired2, json_report, html_report])
        clean_pairs = validate_pairs(read1, read2, args.min_read_pairs)
        # Parse now so a corrupt report cannot silently pass as completed QC.
        json.loads(json_report.read_text())
        qc.update(status="complete", raw_pairs=raw_pairs, clean_pairs=clean_pairs,
                  fastp_json=str(json_report.relative_to(out)), fastp_html=str(html_report.relative_to(out)),
                  unpaired_read1=str(unpaired1.relative_to(out)), unpaired_read2=str(unpaired2.relative_to(out)))
        clean_reads = [read1, read2]
        assembly = work / "assembly" / "spades"
        mode = "--meta" if sample.input_type == "metagenome_reads" else "--isolate"
        command = [*tool("spades"), mode, "-1", str(read1), "-2", str(read2), "-o", str(assembly),
                   "-t", str(threads), "-m", str(memory)]
        run("assembly", [command], [assembly / "contigs.fasta", assembly / "assembly_graph_with_scaffolds.gfa", assembly / "contigs.paths"])
        source = assembly / "contigs.fasta"
        graph = str((assembly / "assembly_graph_with_scaffolds.gfa").relative_to(out))
        contig_paths = str((assembly / "contigs.paths").relative_to(out))
    else:
        source = sample.genome

    canonical = data / "contigs.fna"
    contigs = normalize_fasta(source, canonical, sample.sample_id)
    if clean_reads:
        bam = align(canonical, *clean_reads)
        bai = work / "index" / "reads.bam.bai"
        run("index", [[*tool("samtools"), "index", "-@", str(max(0, threads - 1)), "-o", str(bai), str(bam)]], [bai])
        depth = work / "depth" / "depth.tsv"
        run("depth", [[*tool("depth"), "--outputDepth", str(depth), str(bam)]], [depth])
        coverage = coverage_rows(depth, contigs, sample.sample_id)
    else:
        coverage = [{"sample_id": sample.sample_id, "contig_id": row["contig_id"], "mean_depth": "", "status": "not_available_no_reads"} for row in contigs]

    metagenome = sample.input_type == "metagenome_reads"
    evaluation = work / "assembly_qc" / "report"
    command = [*tool("metaquast" if metagenome else "quast"), str(canonical), "-o", str(evaluation),
               "--threads", str(threads), "--min-contig", "1"]
    if metagenome:
        command += ["--max-ref-number", "0"]
    run("assembly_qc", [command], [evaluation])
    if not any(evaluation.rglob("report.tsv")):
        raise RuntimeError(f"Assembly QC produced no report.tsv for {sample.sample_id}")

    if metagenome and any(row["length"] >= args.min_bin_contig for row in contigs):
        directory = work / "binning"
        prefix = directory / "bin"
        run("binning", [[*tool("metabat2"), "-i", str(canonical), "-a", str(depth), "-o", str(prefix),
                         "-t", str(threads), "-m", str(args.min_bin_contig), "--seed", "42", "--saveCls", "--unbinned"]],
            [directory])
        bins = {f"{sample.sample_id}__bin{int(path.stem.split('.')[-1]):06d}": path
                for path in sorted(directory.glob("bin.*.fa")) if path.stem.split(".")[-1].isdigit()}
    elif metagenome:
        bins = {}
        LOG.warning("[%s] All contigs are below the binning length threshold; preserving all contigs", sample.sample_id)
    else:
        bins = {f"{sample.sample_id}__genome": canonical}
    membership = build_membership(canonical, bins, data / "unbinned.fna", args.min_bin_contig)
    for row in membership:
        row["sample_id"] = sample.sample_id
        if not metagenome:
            row["binning_status"] = "supplied_mag" if sample.input_type == "MAGs" else "isolate_genome"
            row["host_assignment_status"] = "supplied_mag_membership" if sample.input_type == "MAGs" else "isolate_assembly"
    by_contig = {row["contig_id"]: row for row in membership}
    for row in contigs:
        row.update(by_contig[row["contig_id"]])
        row["fasta"] = str(canonical.relative_to(out))
    genomes = [{"sample_id": sample.sample_id, "genome_id": genome_id, "source_fasta": path,
                "genome_type": "mag" if metagenome or sample.input_type == "MAGs" else "isolate"} for genome_id, path in bins.items()]
    if not genomes:
        LOG.warning("[%s] No bins recovered; all %d contigs retained for downstream evidence", sample.sample_id, len(contigs))
    LOG.info("[%s] Prepared %d contigs and %d genomes", sample.sample_id, len(contigs), len(genomes))
    return {"sample": sample, "contigs": contigs, "membership": membership, "genomes": genomes, "coverage": coverage,
            "qc": qc, "assembly_graph": graph, "contig_paths": contig_paths, "assembly_report": str(evaluation.relative_to(out)),
            "reads": [str(path.relative_to(out)) for path in clean_reads],
            "bam": str(bam.relative_to(out)) if bam else "", "bai": str(bai.relative_to(out)) if bai else ""}


def export_results(results, args, runner, tools, databases):
    out = args.outdir
    catalog = out / "catalog"
    genome_dir = catalog / "genomes"
    genome_dir.mkdir(parents=True, exist_ok=True)
    genomes = [row for result in results for row in result["genomes"]]
    for row in genomes:
        destination = genome_dir / (row["genome_id"] + ".fna")
        shutil.copyfile(row.pop("source_fasta"), destination)
        row["fasta"] = str(destination.relative_to(out))
    quality = {}
    env = worker_environment(args.threads)
    if genomes:
        directory = out / "work" / "cohort" / "checkm2"
        report = directory / "report" / "quality_report.tsv"
        runner.run("cohort.checkm2", [[*tools["checkm2"], "predict", "--threads", str(args.threads),
                    "--input", str(genome_dir), "--extension", "fna", "--output-directory", str(report.parent),
                    "--database_path", databases["checkm2"]["path"]]], directory, outputs=[report], env=env)
        quality = parse_checkm2(report)
        expected_ids = {row["genome_id"] for row in genomes}
        if set(quality) != expected_ids:
            raise RuntimeError(f"CheckM2 genome IDs do not match input genomes: missing={sorted(expected_ids - set(quality))}, extra={sorted(set(quality) - expected_ids)}")
    for row in genomes:
        row.update(quality[row["genome_id"]])
        row["comparison_eligible"] = row["completeness"] >= args.min_completeness and row["contamination"] <= args.max_contamination
        row["quality_status"] = "eligible" if row["comparison_eligible"] else "quality_excluded"
    eligible = [row for row in genomes if row["comparison_eligible"]]
    taxonomy = {}
    taxonomy_failures = {}
    if eligible and not args.skip_taxonomy:
        batch = catalog / "gtdbtk.batch.tsv"
        batch.write_text("".join(f"{out / row['fasta']}\t{row['genome_id']}\n" for row in eligible))
        directory = out / "work" / "cohort" / "taxonomy"
        report = directory / "report"
        runner.run("cohort.taxonomy", [[*tools["gtdbtk"], "classify_wf", "--batchfile", str(batch),
                    "--out_dir", str(report), "--cpus", str(args.threads)]],
                   directory, outputs=[report],
                   env={**env, "GTDBTK_DATA_PATH": databases["gtdbtk"]["path"]})
        summaries = sorted(report.glob("gtdbtk.*.summary.tsv"))
        taxonomy = parse_gtdb(summaries)
        taxonomy_failures = collect_gtdb_failures(report)
        expected_ids = {row["genome_id"] for row in eligible}
        extras = (set(taxonomy) | set(taxonomy_failures)) - expected_ids
        if extras:
            raise RuntimeError(f"GTDB-Tk returned unexpected genome IDs: {sorted(extras)}")
        missing = expected_ids - set(taxonomy) - set(taxonomy_failures)
        if missing:
            raise RuntimeError(f"GTDB-Tk omitted genomes without documented failure records: {sorted(missing)}")
    elif not eligible:
        LOG.info("No genomes meet comparison eligibility; retaining all contigs and genome QC results")
    groups = group_genomes(genomes, taxonomy, args.tax_rank)
    for row in groups:
        if args.skip_taxonomy and row["status"] != "quality_excluded":
            row["status"] = "taxonomy_skipped"
        elif row["genome_id"] in taxonomy_failures and row["status"] != "assigned":
            row["status"] = "taxonomy_failed"
    taxonomy_rows = []
    group_status = {row["genome_id"]: row["status"] for row in groups}
    for genome in genomes:
        row = {"sample_id": genome["sample_id"], "genome_id": genome["genome_id"], **taxonomy.get(genome["genome_id"], {})}
        row["status"] = group_status[genome["genome_id"]]
        row["failure_reason"] = taxonomy_failures.get(genome["genome_id"], "")
        taxonomy_rows.append(row)
    contigs = [row for result in results for row in result["contigs"]]
    fields = ["sample_id", "assembly_id", "contig_id", "original_id", "length", "gc_percent", "sequence_sha256", "fasta", "genome_id", "binning_status", "host_assignment_status"]
    write_tsv(out / "contigs.tsv", contigs, fields)
    write_tsv(out / "contig_bin_membership.tsv", [row for result in results for row in result["membership"]],
              ["sample_id", "contig_id", "genome_id", "binning_status", "host_assignment_status"])
    write_tsv(out / "genomes.tsv", genomes, ["sample_id", "genome_id", "genome_type", "fasta", "completeness", "contamination", "comparison_eligible", "quality_status"])
    write_tsv(out / "coverage.tsv", [row for result in results for row in result["coverage"]], ["sample_id", "contig_id", "mean_depth", "status"])
    write_tsv(out / "taxonomy.tsv", taxonomy_rows, ["sample_id", "genome_id", "classification", *RANKS, "classification_method", "warnings", "status", "failure_reason"])
    group_fields = ["sample_id", "genome_id", "tax_rank", "taxon", "group_id", "status"]
    write_tsv(out / "lineage_groups.tsv", groups, group_fields)
    group_dir = out / "groups"
    group_dir.mkdir(exist_ok=True)
    grouped = {}
    for row in groups:
        if row["group_id"]:
            grouped.setdefault(row["group_id"], []).append(row)
    for group_id, members in sorted(grouped.items()):
        write_tsv(group_dir / (group_id + ".tsv"), members, group_fields)
    metadata_fields = sorted({key for result in results for key in result["sample"].metadata})
    sample_rows = [{"sample_id": result["sample"].sample_id, "input_type": result["sample"].input_type,
                    "clean_read1": result["reads"][0] if result["reads"] else "", "clean_read2": result["reads"][1] if result["reads"] else "",
                    "bam": result["bam"], "bai": result["bai"], "assembly_graph": result["assembly_graph"], "contig_paths": result["contig_paths"],
                    "assembly_report": result["assembly_report"], **{f"metadata.{key}": value for key, value in result["sample"].metadata.items()}} for result in results]
    write_tsv(out / "samples.tsv", sample_rows, ["sample_id", "input_type", "clean_read1", "clean_read2", "bam", "bai", "assembly_graph", "contig_paths", "assembly_report", *["metadata." + key for key in metadata_fields]])
    write_tsv(out / "read_qc.tsv", [result["qc"] for result in results], ["sample_id", "status", "raw_pairs", "clean_pairs", "fastp_json", "fastp_html", "unpaired_read1", "unpaired_read2"])
    summary = {"samples": len(results), "contigs": len(contigs), "genomes": len(genomes), "eligible_genomes": len(eligible),
               "unbinned_contigs": sum(not row["genome_id"] for row in contigs), "tax_rank": args.tax_rank,
               "taxonomic_groups": len({row["group_id"] for row in groups if row["group_id"]})}
    write_json(out / "summary.json", summary)
    create_report(out, summary, sample_rows, groups)
    return summary


def create_report(out, summary, samples, groups):
    def table(rows, keys):
        head = "".join(f"<th>{html.escape(key)}</th>" for key in keys)
        body = "".join("<tr>" + "".join(f"<td>{html.escape(str(row.get(key, '')))}</td>" for key in keys) + "</tr>" for row in rows)
        return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    files = ["contigs.tsv", "contig_bin_membership.tsv", "genomes.tsv", "coverage.tsv", "taxonomy.tsv", "lineage_groups.tsv", "read_qc.tsv", "samples.tsv", "run_info.json"]
    links = " ".join(f'<a href="{name}">{name}</a>' for name in files)
    content = '<!doctype html><html lang="en"><meta charset="utf-8"><title>metaMARS preparation report</title>'
    content += '<style>body{font:16px system-ui;margin:3em;max-width:1200px;color:#16313c}table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #ccc;padding:.5em;text-align:left}a{display:inline-block;margin:.4em}th{background:#edf4f5}</style>'
    content += "<h1>metaMARS preparation report</h1>" + table([summary], list(summary))
    content += "<p>All contigs are retained. Bin membership supports host assignment but does not establish a transfer event. Taxonomic groups do not represent abundance or independent biological replicates.</p>"
    content += "<h2>Samples</h2>" + table(samples, ["sample_id", "input_type", "assembly_report"])
    content += "<h2>Genome groups</h2>" + table(groups, ["sample_id", "genome_id", "tax_rank", "taxon", "status"])
    content += "<h2>Outputs</h2>" + links + "</html>"
    (out / "report.html").write_text(content)


def prepare(samples, args):
    """Run independent samples, then evaluate and classify their genomes."""
    out = args.outdir = args.outdir.expanduser().resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError("Output directory must be new or empty. Choose a new --outdir for each run.")
    # Reject mistaken local paths before starting a potentially large download.
    check_path = args.checkm2_db or os.environ.get("CHECKM2DB")
    gtdb_path = args.gtdbtk_db or os.environ.get("GTDBTK_DATA_PATH")
    if check_path:
        validate_checkm2(check_path)
    if gtdb_path and not args.skip_taxonomy:
        validate_gtdbtk(gtdb_path)
    tools = discover_tools(required_tools(samples, args.skip_taxonomy))
    jobs = min(args.jobs, args.threads, args.memory, len(samples))
    threads, memory = args.threads // jobs, args.memory // jobs
    out.mkdir(parents=True, exist_ok=True)
    # Creating work exclusively also prevents two launches sharing an empty directory.
    (out / "work").mkdir()
    handlers = [logging.StreamHandler(), logging.FileHandler(out / "run.log")]
    for handler in handlers:
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S"))
        LOG.addHandler(handler)
    LOG.setLevel(logging.INFO)
    runner = Runner(out, verbose=args.verbose)
    try:
        LOG.info("metaMARS %s: %d samples, %d concurrent jobs, %d threads and %d GB per sample",
                 __version__, len(samples), jobs, threads, memory)
        LOG.info("Taxonomic rank: %s. Memory is advisory outside SPAdes.", args.tax_rank)
        cache = args.db_dir.expanduser().resolve()
        databases = {"checkm2": ensure_checkm2(args.checkm2_db, cache, runner, tools["checkm2"])}
        if not args.skip_taxonomy:
            databases["gtdbtk"] = ensure_gtdbtk(args.gtdbtk_db, cache)
        LOG.info("Database paths: %s", "; ".join(f"{name}={value['path']}" for name, value in databases.items()))
        write_json(out / "run_info.json", {
            "version": __version__, "parameters": asdict(args),
            "samples": [asdict(sample) for sample in samples], "tools": tools, "databases": databases,
        })
        executor = ThreadPoolExecutor(max_workers=jobs)
        try:
            futures = [executor.submit(process_sample, sample, args, runner, tools, threads, memory)
                       for sample in samples]
            results = [future.result() for future in as_completed(futures)]
        except BaseException:
            runner.cancel()
            executor.shutdown(wait=True, cancel_futures=True)
            raise
        else:
            executor.shutdown(wait=True)
        results.sort(key=lambda result: result["sample"].sample_id)
        summary = export_results(results, args, runner, tools, databases)
        LOG.info("Complete: %d genomes, %d contigs, %d taxonomic groups. Report: %s",
                 summary["genomes"], summary["contigs"], summary["taxonomic_groups"], out / "report.html")
    except BaseException as exc:
        runner.cancel()
        LOG.error("Run stopped: %s. Partial outputs and logs are in %s", exc or "interrupted", out)
        raise
    finally:
        for handler in handlers:
            LOG.removeHandler(handler)
            handler.close()
