"""Flat CLI contract and sequencing-format checks; no external tools required."""

from contextlib import redirect_stderr, redirect_stdout
import gzip
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from metamars_prep.cli import input_samples, main, parser
from metamars_prep.pipeline import coverage_rows, fastq_records, validate_pairs
from metamars_prep.taxonomy import RANKS


class CliTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.genome = self.root / "genome.fa"
        self.genome.write_text(">contig\nACGT\n")
        self.outdir = self.root / "results"
        self.sheet = self.root / "samples.tsv"
        self.sheet.write_text("sample_id\tinput_type\tgenome\nsample1\tisolate_genomes\tgenome.fa\n")
        self.single = ["--samples", str(self.sheet), "--outdir", str(self.outdir)]

    def test_public_flags_are_exactly_the_requested_flat_interface(self):
        flags = {flag for action in parser()._actions for flag in action.option_strings}
        self.assertEqual(flags, {
            "-h", "--help", "-v", "--verbose", "-o", "--outdir", "-t", "--threads", "-j", "--jobs",
            "-m", "--memory", "-tr", "--tax-rank", "-st", "--skip-taxonomy",
            "--db-dir", "--checkm2-db", "--gtdbtk-db", "-it", "--input-type", "-s", "--samples",
        })
        args = parser().parse_args(self.single)
        self.assertEqual((args.threads, args.jobs, args.memory), (2, 1, 8))
        self.assertEqual(args.tax_rank, "species")
        self.assertFalse(args.skip_taxonomy)
        self.assertFalse(args.verbose)
        self.assertTrue(parser().parse_args(self.single + ["-v"]).verbose)
        self.assertTrue(parser().parse_args(self.single + ["--verbose"]).verbose)
        self.assertIsNone(args.input_type)
        self.assertEqual(args.db_dir, Path.home() / "databases" / "metamars")
        self.assertIsNone(args.checkm2_db)
        self.assertIsNone(args.gtdbtk_db)
        self.assertEqual(args.samples, self.sheet)
        self.assertFalse(hasattr(args, "command"))

    def test_help_shows_defaults_and_succeeds_without_required_inputs(self):
        output = io.StringIO()
        with redirect_stdout(output), self.assertRaises(SystemExit) as exit_status:
            main(["-h"])
        self.assertEqual(exit_status.exception.code, 0)
        help_text = output.getvalue()
        for default in ("dflt=2", "dflt=1", "dflt=8", "dflt=species", "dflt=False"):
            self.assertIn(default, help_text)
        self.assertIn("--gtdbtk-db", help_text)
        self.assertNotIn("--dry-run", help_text)
        self.assertFalse(self.outdir.exists())

    def test_short_flags_remain_distinct_from_their_prefixes(self):
        args = parser().parse_args([
            "-s", str(self.sheet), "-o", str(self.outdir), "-t", "4", "-tr", "genus",
            "-j", "2", "-m", "16", "-st", "-it", "isolate_genomes",
            "--gtdbtk-db", str(self.root / "gtdb"),
        ])
        self.assertEqual((args.threads, args.jobs, args.memory), (4, 2, 16))
        self.assertEqual(args.tax_rank, "genus")
        self.assertEqual(args.samples, self.sheet)
        self.assertTrue(args.skip_taxonomy)
        self.assertEqual(args.input_type, "isolate_genomes")

    def test_samples_and_outdir_are_both_required(self):
        for arguments in ([], ["-s", str(self.sheet)], ["-o", str(self.outdir)]):
            with self.subTest(arguments=arguments), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    parser().parse_args(arguments)
                self.assertEqual(error.exception.code, 2)

    def test_all_taxonomic_ranks_supported_without_lineage_flags(self):
        for rank in RANKS:
            with self.subTest(rank=rank):
                self.assertEqual(parser().parse_args(self.single + ["--tax-rank", rank]).tax_rank, rank)
        for options in [["--tax-rank", "strain"], ["--lineage", "Escherichia"],
                        ["--lineage-rank", "genus"], ["--lineages", "Escherichia"]]:
            with self.subTest(options=options), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    parser().parse_args(self.single + options)
                self.assertEqual(error.exception.code, 2)

    def test_removed_flags_and_subcommands_are_rejected(self):
        cases = [[flag] for flag in ("prepare", "doctor", "--version", "--dry-run", "--quiet", "--resume", "-r")]
        cases += [[flag, "value"] for flag in (
            "--memory-gb", "--sample-id", "--read1", "--read2", "--genome",
            "--tool", "--host-reference", "--reference", "--min-bin-contig", "--min-completeness",
            "--max-contamination", "--min-read-length", "--qualified-phred", "--min-read-pairs",
            "--adapter-fasta", "--tax", "--sam",
        )]
        cases += [["--allow-reference-download"]]
        for options in cases:
            with self.subTest(options=options), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    parser().parse_args(self.single + options)
                self.assertEqual(error.exception.code, 2)

    def test_numeric_options_reject_invalid_values(self):
        cases = [("--threads", "0"), ("--threads", "word"), ("--jobs", "-1"),
                 ("--memory", "0"), ("--memory", "2.5"), ("--memory", "8GB")]
        for option, value in cases:
            with self.subTest(option=option, value=value), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    parser().parse_args(self.single + [option, value])
                self.assertEqual(error.exception.code, 2)

    def test_input_type_choices_are_exact(self):
        for input_type in ("metagenome_reads", "isolate_reads", "isolate_genomes", "MAGs"):
            self.assertEqual(parser().parse_args(self.single + ["-it", input_type]).input_type, input_type)
        for input_type in ("isolate_genome", "mag", "mags", "MAG", "reads"):
            with self.subTest(input_type=input_type), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    parser().parse_args(self.single + ["-it", input_type])

    def test_default_input_type_applies_to_absent_or_blank_sheet_values(self):
        for contents in (
            "sample_id\tgenome\nsample1\tgenome.fa\n",
            "sample_id\tinput_type\tgenome\nsample1\t\tgenome.fa\n",
        ):
            with self.subTest(contents=contents):
                self.sheet.write_text(contents)
                samples = input_samples(parser().parse_args(self.single + ["-it", "isolate_genomes"]))
                self.assertEqual(len(samples), 1)
                self.assertEqual(samples[0].input_type, "isolate_genomes")
                self.assertEqual(samples[0].genome, self.genome)

    def test_conflicting_sheet_type_fails_before_preflight(self):
        with patch("metamars_prep.pipeline.discover_tools") as discovery, redirect_stderr(io.StringIO()) as error:
            result = main(self.single + ["-it", "MAGs"])
        self.assertEqual(result, 1)
        self.assertIn("input_type", error.getvalue())
        discovery.assert_not_called()
        self.assertFalse(self.outdir.exists())

    def test_missing_sheet_type_without_default_creates_nothing(self):
        self.sheet.write_text("sample_id\tgenome\nsample1\tgenome.fa\n")
        with patch("metamars_prep.pipeline.discover_tools") as discovery, redirect_stderr(io.StringIO()) as error:
            result = main(self.single)
        self.assertEqual(result, 1)
        self.assertIn("input_type", error.getvalue())
        discovery.assert_not_called()
        self.assertFalse(self.outdir.exists())

    def test_missing_tool_fails_before_output_creation(self):
        with patch("metamars_prep.pipeline.discover_tools", side_effect=ValueError("Missing tools: checkm2")), redirect_stderr(io.StringIO()) as error:
            result = main(self.single)
        self.assertEqual(result, 1)
        self.assertIn("Missing tools", error.getvalue())
        self.assertFalse(self.outdir.exists())

    def test_gzip_eof_is_reported_without_traceback(self):
        with patch("metamars_prep.pipeline.prepare", side_effect=EOFError("Compressed file ended early")), redirect_stderr(io.StringIO()) as error:
            self.assertEqual(main(self.single), 1)
        self.assertIn("Compressed file ended early", error.getvalue())
        self.assertNotIn("Traceback", error.getvalue())

    def test_external_tool_output_path_restrictions_fail_before_writing(self):
        for name in ("output with spaces", "output;echo", "output$(echo)"):
            with self.subTest(name=name), redirect_stderr(io.StringIO()) as error:
                destination = self.root / name
                result = main(self.single + ["--outdir", str(destination)])
                self.assertEqual(result, 1)
                self.assertIn("external tools", error.getvalue())
                self.assertFalse(destination.exists())


class SequenceHelperTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def reads(self, name, content):
        path = self.root / name
        if path.suffix == ".gz":
            with gzip.open(path, "wt") as stream:
                stream.write(content)
        else:
            path.write_text(content)
        return path

    def test_paired_gzip_and_illumina_headers_are_supported(self):
        read1 = self.reads("R1.fq.gz", "@first/1\nACGT\n+\nIIII\n@second 1:N:0:1\nNNNN\n+\n!!!!\n")
        read2 = self.reads("R2.fq.gz", "@first/2\nACGT\n+\nIIII\n@second 2:N:0:1\nNNNN\n+\n!!!!\n")
        self.assertEqual(validate_pairs(read1, read2, 2), 2)

    def test_mismatched_ids_counts_and_insufficient_pairs_fail(self):
        read1 = self.reads("R1.fq", "@first/1\nACGT\n+\nIIII\n")
        for content in ["@other/2\nACGT\n+\nIIII\n", "", "@first/2\nACGT\n+\nIIII\n@second/2\nACGT\n+\nIIII\n"]:
            with self.subTest(content=content):
                read2 = self.reads("R2.fq", content)
                with self.assertRaisesRegex(ValueError, "Mismatched"):
                    validate_pairs(read1, read2, 1)
        read2 = self.reads("R2.fq", "@first/2\nACGT\n+\nIIII\n")
        with self.assertRaisesRegex(ValueError, "required at least 2"):
            validate_pairs(read1, read2, 2)

    def test_malformed_truncated_and_empty_name_records_fail(self):
        for content in ["@first\nACGT\n+\n", "first\nACGT\n+\nIIII\n",
                        "@first\nACGU\n+\nIIII\n", "@first\nACGT\n+\nII I\n",
                        "@\nACGT\n+\nIIII\n", "@first\nACGT\nnotplus\nIIII\n"]:
            with self.subTest(content=content):
                path = self.reads("invalid.fq.gz", content)
                with self.assertRaises(ValueError):
                    list(fastq_records(path))

    def test_coverage_distinguishes_missing_measurement_and_observed_zero(self):
        path = self.root / "depth.tsv"
        path.write_text("contigName\ttotalAvgDepth\na\t0\n")
        rows = coverage_rows(path, [{"contig_id": "a"}, {"contig_id": "b"}], "sample1")
        self.assertEqual((rows[0]["mean_depth"], rows[0]["status"]), (0, "measured"))
        self.assertEqual((rows[1]["mean_depth"], rows[1]["status"]), ("", "unavailable"))
        for body in ["a\tnan\n", "a\t-1\n", "a\t1\na\t1\n", "unknown\t1\n", "a\n", "a\t\n"]:
            with self.subTest(body=body):
                path.write_text("contigName\ttotalAvgDepth\n" + body)
                with self.assertRaises(ValueError):
                    coverage_rows(path, [{"contig_id": "a"}], "sample1")


if __name__ == "__main__":
    unittest.main()
