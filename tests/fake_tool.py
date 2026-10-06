"""Synthetic external tools for CLI integration tests.

These fixtures validate orchestration, file contracts, and provenance only.
They do not model bioinformatics algorithms or validate scientific accuracy.
The tests install copies under tool names with the current Python shebang.
"""

import gzip
import json
import os
from pathlib import Path
import sys


def option(arguments, name):
    return arguments[arguments.index(name) + 1]


def write(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def read_fasta(path):
    records = []
    name, sequence = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(sequence)))
            name, sequence = line[1:].split()[0], []
        else:
            sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence)))
    return records


def fasta_text(records):
    return "".join(f">{name}\n{sequence}\n" for name, sequence in records)


def read_fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return handle.read()


def write_fastq(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        handle.write(text)


def main():
    name = Path(sys.argv[0]).name
    name = {"spades.py": "spades", "quast.py": "quast", "metaquast.py": "metaquast", "jgi_summarize_bam_contig_depths": "depth"}.get(name, name)
    arguments = sys.argv[1:]
    if arguments in (["--version"], ["--help"]):
        print(f"{name}: synthetic integration fixture 1.0")
        return 0
    call_log = os.environ.get("METAMARS_FAKE_CALL_LOG")
    if call_log:
        record = {"tool": name, "arguments": arguments, "cwd": str(Path.cwd())}
        descriptor = os.open(call_log, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        try:
            os.write(descriptor, (json.dumps(record) + "\n").encode())
        finally:
            os.close(descriptor)
    print(f"Synthetic {name} fixture started", file=sys.stderr, flush=True)
    if os.environ.get("METAMARS_FAKE_FAIL") == name:
        print(f"Requested synthetic failure in {name}", file=sys.stderr, flush=True)
        return 17
    if name == "fastp":
        for source, destination in (("--in1", "--out1"), ("--in2", "--out2")):
            write_fastq(option(arguments, destination), read_fastq(option(arguments, source)))
        for unpaired in ("--unpaired1", "--unpaired2"):
            write_fastq(option(arguments, unpaired), "")
        write(option(arguments, "--json"), json.dumps({"summary": {"before_filtering": {"total_reads": 2}, "after_filtering": {"total_reads": 2}}}))
        write(option(arguments, "--html"), "<!doctype html><title>Synthetic fastp report</title>")
    elif name == "spades":
        directory = Path(option(arguments, "-o"))
        repeat = 250 if os.environ.get("METAMARS_FAKE_SHORT") == "1" else 750
        write(directory / "contigs.fasta", fasta_text([("NODE_1_length_3000_cov_20", "ACGT" * repeat), ("NODE_2_length_600_cov_5", "GATTAC" * 100)]))
        write(directory / "assembly_graph_with_scaffolds.gfa", "H\tVN:Z:1.0\nS\t1\tACGT\n")
        write(directory / "contigs.paths", "NODE_1_length_3000_cov_20\n1+\nNODE_2_length_600_cov_5\n2+\n")
    elif name == "minimap2":
        records = read_fasta(arguments[-3])
        sam = "@HD\tVN:1.6\tSO:unsorted\n" + "".join(f"@SQ\tSN:{identifier}\tLN:{len(sequence)}\n" for identifier, sequence in records)
        if "-o" in arguments:
            write(option(arguments, "-o"), sam)
        else:
            print(sam, end="", flush=True)
    elif name == "samtools":
        if arguments[0] == "sort":
            content = sys.stdin.read() if arguments[-1] == "-" else Path(arguments[-1]).read_text()
            write(option(arguments, "-o"), content)
        elif arguments[0] == "index":
            write(option(arguments, "-o"), "synthetic index\n")
        elif arguments[0] == "fastq":
            for mate in (1, 2):
                write_fastq(option(arguments, f"-{mate}"), f"@pair1/{mate}\n{'ACGT' * 25}\n+\n{'I' * 100}\n")
        else:
            raise ValueError(f"Unsupported synthetic samtools operation: {arguments}")
    elif name == "depth":
        rows = ["contigName\tcontigLen\ttotalAvgDepth\treads.bam\treads.bam-var"]
        for line in Path(arguments[-1]).read_text().splitlines():
            if line.startswith("@SQ\t"):
                fields = dict(item.split(":", 1) for item in line.split("\t")[1:])
                rows.append(f"{fields['SN']}\t{fields['LN']}\t20\t20\t0")
        write(option(arguments, "--outputDepth"), "\n".join(rows) + "\n")
    elif name in {"quast", "metaquast"}:
        directory = Path(option(arguments, "-o"))
        write(directory / "report.tsv", "Assembly\tsynthetic\n# contigs\t2\n")
        write(directory / "report.html", "<!doctype html><title>Synthetic assembly report</title>")
    elif name == "metabat2":
        if os.environ.get("METAMARS_FAKE_ZERO_BINS") != "1" and "nobins" not in str(Path.cwd()):
            records = read_fasta(option(arguments, "-i"))
            eligible = [record for record in records if len(record[1]) >= int(option(arguments, "-m"))]
            if eligible:
                write(option(arguments, "-o") + ".1.fa", fasta_text(eligible[:1]))
    elif name == "checkm2":
        directory = Path(option(arguments, "--input"))
        rows = ["Name\tCompleteness\tContamination"]
        for genome in sorted(directory.glob("*.fna")):
            completeness, contamination = (25, 20) if "quality_bad" in genome.stem else (98, 1)
            rows.append(f"{genome.stem}\t{completeness}\t{contamination}")
        write(Path(option(arguments, "--output-directory")) / "quality_report.tsv", "\n".join(rows) + "\n")
    elif name == "gtdbtk":
        directory = Path(option(arguments, "--out_dir"))
        directory.mkdir(parents=True, exist_ok=True)
        if os.environ.get("METAMARS_FAKE_GTDB") == "empty":
            return 0
        if os.environ.get("METAMARS_FAKE_GTDB") == "failed":
            failed = []
            for line in Path(option(arguments, "--batchfile")).read_text().splitlines():
                _, genome = line.split("\t")
                failed.append(f"{genome}\tSynthetic insufficient marker genes")
            write(directory / "gtdbtk.failed_genomes.tsv", "\n".join(failed) + "\n")
            return 0
        rows = ["user_genome\tclassification\tclassification_method\twarnings"]
        for line in Path(option(arguments, "--batchfile")).read_text().splitlines():
            _, genome = line.split("\t")
            if os.environ.get("METAMARS_FAKE_GTDB") == "missing":
                continue
            species = "Escherichia fergusonii" if "iso_b" in genome else "Escherichia coli"
            taxonomy = f"d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__{species}"
            rows.append(f"{genome}\t{taxonomy}\tsynthetic_fixture\t")
        write(directory / "gtdbtk.bac120.summary.tsv", "\n".join(rows) + "\n")
    else:
        raise ValueError(f"Unknown synthetic tool name: {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
