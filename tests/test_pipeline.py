"""CLI integration tests with explicitly synthetic external-tool fixtures.

These tests run Python orchestration and output contracts end to end. They
are not scientific validation of SPAdes, MetaBAT2, CheckM2, QUAST or GTDB-Tk.
No external tools, reference downloads, or network access are required.
"""

from collections import Counter
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


PROJECT = Path(__file__).resolve().parents[1]
TOOLS = ("fastp", "spades.py", "minimap2", "samtools", "jgi_summarize_bam_contig_depths", "metabat2", "metaquast.py", "quast.py", "checkm2", "gtdbtk")


class PipelineIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="metamars-integration-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.out = self.root / "results"
        self.bin = self.root / "fake_bin"
        self.bin.mkdir()
        body = Path(__file__).with_name("fake_tool.py").read_text()
        for tool in TOOLS:
            executable = self.bin / tool
            executable.write_text(f"#!{sys.executable}\n" + body)
            executable.chmod(0o755)
        self.call_log = self.root / "tool_calls.jsonl"
        self.environment = {**os.environ, "METAMARS_FAKE_CALL_LOG": str(self.call_log), "PYTHONDONTWRITEBYTECODE": "1"}
        self.environment.pop("METAMARS_FAKE_FAIL", None)
        self.environment.pop("METAMARS_FAKE_ZERO_BINS", None)
        self.environment.pop("METAMARS_FAKE_GTDB", None)
        self.environment.pop("METAMARS_FAKE_SHORT", None)
        self.environment["PATH"] = str(self.bin) + os.pathsep + os.environ.get("PATH", "")
        self.checkm_db = self.root / "checkm2.dmnd"
        self.checkm_db.write_text("synthetic database fixture\n")
        self.environment["CHECKM2DB"] = str(self.checkm_db)
        self.gtdb_db = self.root / "gtdb"
        self.gtdb_db.mkdir()
        (self.gtdb_db / "VERSION_DATA").write_text("synthetic database fixture\n")
        self.genome = self.root / "isolate.fna"
        self.genome.write_text(">isolate_contig\n" + "ACGT" * 750 + "\n")
        self.reads = []
        for mate in (1, 2):
            path = self.root / f"reads.R{mate}.fastq.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(f"@pair1/{mate}\n{'ACGT' * 25}\n+\n{'I' * 100}\n")
            self.reads.append(path)

    def invoke(self, arguments, expected=0, env=None, dependencies=True, threads=2):
        command = [sys.executable, "-m", "metamars_prep", "--outdir", str(self.out), "--threads", str(threads)]
        if dependencies:
            command += ["--gtdbtk-db", str(self.gtdb_db)]
        result = subprocess.run(command + arguments, cwd=PROJECT, env={**self.environment, **(env or {})}, text=True, capture_output=True, timeout=60)
        self.assertEqual(result.returncode, expected, msg=f"Command: {command + arguments}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")
        return result

    def single_genome(self, identifier="iso_a", input_type="isolate_genomes"):
        return self.sheet([{"sample_id": identifier, "input_type": input_type, "genome": str(self.genome)}])

    def single_reads(self, identifier="meta_a", input_type="metagenome_reads"):
        return self.sheet([{"sample_id": identifier, "input_type": input_type, "read1": str(self.reads[0]), "read2": str(self.reads[1])}])

    def rows(self, filename):
        with (self.out / filename).open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def calls(self):
        if not self.call_log.exists():
            return []
        return [json.loads(line) for line in self.call_log.read_text().splitlines()]

    def info(self):
        return json.loads((self.out / "run_info.json").read_text())

    def summary(self):
        return json.loads((self.out / "summary.json").read_text())

    def sheet(self, rows):
        sheet = self.root / "samples.tsv"
        fields = ["sample_id", "input_type", "read1", "read2", "genome", "subject_id", "cohort"]
        with sheet.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return ["--samples", str(sheet)]

    def test_single_isolate_genome_outputs_quality_taxonomy_and_run_info(self):
        result = self.invoke(self.single_genome() + ["--tax-rank", "genus"])
        genomes = self.rows("genomes.tsv")
        self.assertEqual(len(genomes), 1)
        self.assertEqual(genomes[0]["genome_id"], "iso_a__genome")
        self.assertEqual(float(genomes[0]["completeness"]), 98)
        self.assertEqual(float(genomes[0]["contamination"]), 1)
        groups = self.rows("lineage_groups.tsv")
        self.assertEqual(groups[0]["tax_rank"], "genus")
        self.assertEqual(groups[0]["taxon"], "g__Escherichia")
        self.assertTrue(groups[0]["group_id"])
        self.assertEqual({call["tool"] for call in self.calls()}, {"quast", "checkm2", "gtdbtk"})
        info = self.info()
        self.assertEqual(info["tools"]["quast"], [str(self.bin / "quast.py")])
        self.assertEqual(info["databases"]["gtdbtk"]["path"], str(self.gtdb_db))
        self.assertFalse((self.out / ".metamars").exists())
        self.assertFalse((self.out / "run_manifest.json").exists())
        self.assertTrue((self.out / "report.html").is_file())
        self.assertFalse((self.out / ".metamars.lock").exists())
        self.assertIn("Synthetic quast fixture started", result.stderr)

    def test_short_options_and_global_type_reach_workflow_resources(self):
        sheet = self.root / "homogeneous.csv"
        sheet.write_text(f"sample_id,read1,read2\nshort_flags,{self.reads[0].name},{self.reads[1].name}\n")
        self.invoke(["-s", str(sheet), "-it", "metagenome_reads", "-t", "4", "-j", "1", "-m", "12", "-tr", "genus", "-st"])
        assembly = next(call["arguments"] for call in self.calls() if call["tool"] == "spades")
        self.assertEqual(assembly[assembly.index("-t") + 1], "4")
        self.assertEqual(assembly[assembly.index("-m") + 1], "12")
        self.assertNotIn("gtdbtk", {call["tool"] for call in self.calls()})
        self.assertEqual(self.info()["parameters"]["memory"], 12)
        self.assertEqual(self.rows("samples.tsv")[0]["input_type"], "metagenome_reads")
        self.assertEqual(self.rows("lineage_groups.tsv")[0]["tax_rank"], "genus")

    def test_metagenome_retains_binned_and_short_unbinned_contigs(self):
        self.invoke(self.single_reads() + ["--tax-rank", "genus"])
        contigs = self.rows("contigs.tsv")
        self.assertEqual(len(contigs), 2)
        by_length = {int(row["length"]): row for row in contigs}
        self.assertEqual(by_length[3000]["genome_id"], "meta_a__bin000001")
        self.assertEqual(by_length[600]["genome_id"], "")
        unbinned = (self.out / "samples/meta_a/unbinned.fna").read_text()
        self.assertIn(by_length[600]["contig_id"], unbinned)
        self.assertNotIn(by_length[3000]["contig_id"], unbinned)
        coverage = self.rows("coverage.tsv")
        self.assertEqual({row["contig_id"] for row in coverage}, {row["contig_id"] for row in contigs})
        self.assertTrue(all(float(row["mean_depth"]) == 20 for row in coverage))
        sample = self.rows("samples.tsv")[0]
        for field in ("bam", "bai", "assembly_graph", "contig_paths", "clean_read1", "clean_read2"):
            self.assertTrue((self.out / sample[field]).is_file(), field)
        calls = self.calls()
        self.assertIn("--meta", next(call["arguments"] for call in calls if call["tool"] == "spades"))
        self.assertIn("--max-ref-number", next(call["arguments"] for call in calls if call["tool"] == "metaquast"))
        self.assertEqual(self.summary()["unbinned_contigs"], 1)

    def test_zero_bins_is_success_and_preserves_all_contigs(self):
        self.invoke(self.single_reads("nobins"))
        self.assertEqual(self.rows("genomes.tsv"), [])
        self.assertEqual(self.rows("lineage_groups.tsv"), [])
        self.assertEqual(len(self.rows("contigs.tsv")), 2)
        self.assertTrue(all(not row["genome_id"] for row in self.rows("contigs.tsv")))
        self.assertNotIn("checkm2", {call["tool"] for call in self.calls()})
        self.assertNotIn("gtdbtk", {call["tool"] for call in self.calls()})
        self.assertTrue((self.out / "summary.json").is_file())

    def test_all_contigs_below_binning_threshold_skip_binner_and_are_retained(self):
        self.invoke(self.single_reads(), env={"METAMARS_FAKE_SHORT": "1"})
        self.assertEqual(self.rows("genomes.tsv"), [])
        self.assertEqual(len(self.rows("contigs.tsv")), 2)
        self.assertTrue(all(row["binning_status"] == "too_short" for row in self.rows("contigs.tsv")))
        self.assertNotIn("metabat2", {call["tool"] for call in self.calls()})

    def test_two_genomes_group_by_genus_or_species_in_separate_runs(self):
        arguments = self.sheet([
            {"sample_id": "iso_a", "input_type": "isolate_genomes", "genome": self.genome.name, "subject_id": "subject_a", "cohort": "case"},
            {"sample_id": "iso_b", "input_type": "isolate_genomes", "genome": self.genome.name, "subject_id": "subject_b", "cohort": "control"},
        ])
        self.invoke(arguments + ["--tax-rank", "genus", "--jobs", "2"])
        self.assertEqual(len({row["group_id"] for row in self.rows("lineage_groups.tsv")}), 1)
        self.assertEqual(len(list((self.out / "groups").glob("*.tsv"))), 1)
        self.out = self.root / "species_results"
        self.invoke(arguments + ["--tax-rank", "species", "--jobs", "2"])
        self.assertEqual({row["taxon"] for row in self.rows("lineage_groups.tsv")}, {"s__Escherichia coli", "s__Escherichia fergusonii"})
        self.assertEqual(len(list((self.out / "groups").glob("*.tsv"))), 2)
        self.assertEqual(self.summary()["tax_rank"], "species")

    def test_mixed_sample_sheet_preserves_metadata_and_uses_distinct_routes(self):
        arguments = self.sheet([
            {"sample_id": "meta_a", "input_type": "metagenome_reads", "read1": self.reads[0].name, "read2": self.reads[1].name, "subject_id": "participant_1", "cohort": "case"},
            {"sample_id": "iso_b", "input_type": "isolate_genomes", "genome": self.genome.name, "subject_id": "participant_2", "cohort": "control"},
            {"sample_id": "mag_a", "input_type": "MAGs", "genome": self.genome.name, "subject_id": "participant_3", "cohort": "case"},
        ])
        self.invoke(arguments + ["--jobs", "2", "--tax-rank", "genus"])
        self.assertEqual(len(self.rows("genomes.tsv")), 3)
        self.assertEqual(len(self.rows("contigs.tsv")), 4)
        samples = {row["sample_id"]: row for row in self.rows("samples.tsv")}
        self.assertEqual(samples["meta_a"]["metadata.subject_id"], "participant_1")
        self.assertEqual(samples["iso_b"]["metadata.cohort"], "control")
        self.assertEqual(Counter(call["tool"] for call in self.calls())["metabat2"], 1)
        self.assertEqual(Counter(call["tool"] for call in self.calls())["quast"], 2)

    def test_isolate_reads_use_isolate_mode_and_single_thread_mapping(self):
        self.invoke(self.single_reads("isolate_reads", "isolate_reads"), threads=1)
        calls = self.calls()
        spades = next(call for call in calls if call["tool"] == "spades")
        self.assertIn("--isolate", spades["arguments"])
        self.assertNotIn("--meta", spades["arguments"])
        self.assertNotIn("metabat2", {call["tool"] for call in calls})
        mapping = next(call for call in calls if call["tool"] == "minimap2")
        self.assertNotIn("-o", mapping["arguments"])
        self.assertEqual(mapping["arguments"][mapping["arguments"].index("-t") + 1], "1")
        self.assertEqual(len(self.rows("genomes.tsv")), 1)

    def test_skip_taxonomy_retains_genome_with_explicit_unassigned_status(self):
        self.invoke(self.single_genome() + ["--skip-taxonomy", "--tax-rank", "genus"])
        self.assertNotIn("gtdbtk", {call["tool"] for call in self.calls()})
        self.assertNotIn("gtdbtk", self.info()["databases"])
        group = self.rows("lineage_groups.tsv")[0]
        self.assertEqual(group["status"], "taxonomy_skipped")
        self.assertEqual(group["group_id"], "")
        self.assertEqual(self.rows("taxonomy.tsv")[0]["status"], "taxonomy_skipped")
        self.assertEqual(len(self.rows("contigs.tsv")), 1)
        self.assertEqual(len(self.rows("genomes.tsv")), 1)

    def test_quality_excluded_genome_keeps_contigs_and_skips_taxonomy(self):
        self.invoke(self.single_genome("quality_bad"))
        genome = self.rows("genomes.tsv")[0]
        self.assertEqual(genome["quality_status"], "quality_excluded")
        self.assertEqual(float(genome["completeness"]), 25)
        self.assertEqual(float(genome["contamination"]), 20)
        self.assertNotIn("gtdbtk", {call["tool"] for call in self.calls()})
        self.assertEqual(self.rows("lineage_groups.tsv")[0]["status"], "quality_excluded")
        self.assertEqual(self.rows("lineage_groups.tsv")[0]["group_id"], "")
        self.assertEqual(self.rows("contigs.tsv")[0]["genome_id"], genome["genome_id"])
        self.assertTrue((self.out / genome["fasta"]).is_file())

    def test_documented_gtdb_failure_retains_genome_and_reason(self):
        self.invoke(self.single_genome(), env={"METAMARS_FAKE_GTDB": "failed"})
        taxonomy = self.rows("taxonomy.tsv")[0]
        self.assertEqual(taxonomy["status"], "taxonomy_failed")
        self.assertEqual(taxonomy["failure_reason"], "Synthetic insufficient marker genes")
        group = self.rows("lineage_groups.tsv")[0]
        self.assertEqual(group["status"], "taxonomy_failed")
        self.assertEqual(group["group_id"], "")
        self.assertEqual(len(self.rows("contigs.tsv")), 1)
        self.assertEqual(len(self.rows("genomes.tsv")), 1)
        self.assertTrue((self.out / "summary.json").is_file())

    def test_tool_failure_keeps_logs_and_requires_new_output_directory(self):
        arguments = self.single_reads()
        result = self.invoke(arguments, expected=1, env={"METAMARS_FAKE_FAIL": "spades"})
        self.assertIn("17", result.stderr)
        self.assertIn("spades", result.stderr)
        self.assertFalse((self.out / "summary.json").exists())
        self.assertFalse((self.out / ".metamars").exists())
        self.assertIn("Run stopped", (self.out / "run.log").read_text())
        calls_before = self.calls()
        rerun = self.invoke(arguments, expected=1)
        self.assertIn("new or empty", rerun.stderr)
        self.assertEqual(self.calls(), calls_before)
        self.out = self.root / "retry"
        self.invoke(arguments)
        counts = Counter(call["tool"] for call in self.calls())
        self.assertEqual(counts["fastp"], 2)
        self.assertEqual(counts["spades"], 2)

    def test_existing_output_is_preserved_without_running_tools(self):
        self.out.mkdir()
        marker = self.out / "keep.txt"
        marker.write_text("user content")
        self.invoke(self.single_genome(), expected=1)
        self.assertEqual(marker.read_text(), "user content")
        self.assertEqual(self.calls(), [])

    def test_empty_existing_output_is_accepted(self):
        self.out.mkdir()
        self.invoke(self.single_genome())
        self.assertEqual(self.summary()["genomes"], 1)

    def test_missing_taxonomy_results_fail_instead_of_silent_success(self):
        for mode in ("empty", "missing"):
            with self.subTest(mode=mode):
                self.out = self.root / mode
                self.invoke(self.single_genome(), expected=1, env={"METAMARS_FAKE_GTDB": mode})
                self.assertFalse((self.out / "summary.json").exists())
                self.assertFalse((self.out / ".metamars.lock").exists())

    def test_deprecated_lineage_options_are_rejected(self):
        result = self.invoke(self.single_genome() + ["--lineage-rank", "genus"], expected=2)
        self.assertIn("unrecognized arguments", result.stderr)
        self.assertFalse(self.out.exists())
        result = self.invoke(self.single_genome() + ["--lineage", "g__Escherichia"], expected=2)
        self.assertIn("unrecognized arguments", result.stderr)


if __name__ == "__main__":
    unittest.main()
