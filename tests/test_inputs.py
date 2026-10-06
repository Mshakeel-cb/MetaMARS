import tempfile
import unittest
from pathlib import Path

from metamars_prep.models import Sample, load_samples, validate_sample


class InputTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.r1 = self.root / "one.fastq"
        self.r2 = self.root / "two.fastq"
        self.genome = self.root / "genome.fna"
        self.r1.write_text("@read/1\nACGT\n+\nIIII\n")
        self.r2.write_text("@read/2\nACGT\n+\nIIII\n")
        self.genome.write_text(">contig\nACGT\n")

    def manifest(self, content, name="samples.tsv"):
        path = self.root / name
        path.write_text(content)
        return path

    def test_mixed_inputs_resolve_relative_paths_and_preserve_metadata(self):
        path = self.manifest(
            "sample_id\tinput_type\tread1\tread2\tgenome\tsubject_id\tgroup\n"
            "m1\tmetagenome_reads\tone.fastq\ttwo.fastq\t\tpatient1\tcontrol\n"
            "i1\tisolate_reads\tone.fastq\ttwo.fastq\t\tpatient2\tcase\n"
            "i2\tisolate_genomes\t\t\tgenome.fna\tpatient3\tcase\n"
            "g1\tMAGs\t\t\tgenome.fna\t\tcontrol\n"
        )
        samples = load_samples(path)
        self.assertEqual(len(samples), 4)
        self.assertEqual(samples[0].files, [self.r1, self.r2])
        self.assertEqual(samples[2].files, [self.genome])
        self.assertEqual(samples[0].metadata, {"subject_id": "patient1", "group": "control"})
        self.assertEqual([sample.input_type for sample in samples], ["metagenome_reads", "isolate_reads", "isolate_genomes", "MAGs"])

    def test_csv_handles_quoted_metadata(self):
        path = self.manifest('sample_id,input_type,genome,location\nisolate,isolate_genomes,genome.fna,"city, region"\n', "samples.csv")
        self.assertEqual(load_samples(path)[0].metadata["location"], "city, region")

    def test_global_type_supplies_omitted_column_and_preserves_metadata(self):
        path = self.manifest("sample_id\tgenome\tgroup\ni1\tgenome.fna\tcase\n")
        sample = load_samples(path, "isolate_genomes")[0]
        self.assertEqual(sample.input_type, "isolate_genomes")
        self.assertEqual(sample.genome, self.genome)
        self.assertEqual(sample.metadata, {"group": "case"})
        with self.assertRaisesRegex(ValueError, "requires an input_type column"):
            load_samples(path)

    def test_global_type_fills_blank_cells_and_accepts_matching_cells(self):
        path = self.manifest(
            "sample_id\tinput_type\tread1\tread2\n"
            "s1\t\tone.fastq\ttwo.fastq\n"
            "s2\tmetagenome_reads\tone.fastq\ttwo.fastq\n"
        )
        samples = load_samples(path, "metagenome_reads")
        self.assertEqual([sample.input_type for sample in samples], ["metagenome_reads", "metagenome_reads"])
        self.assertEqual(samples[0].files, [self.r1, self.r2])
        with self.assertRaisesRegex(ValueError, "input_type is required"):
            load_samples(path)

    def test_global_type_rejects_conflicting_rows(self):
        path = self.manifest("sample_id\tinput_type\tgenome\ns\tMAGs\tgenome.fna\n")
        with self.assertRaisesRegex(ValueError, "conflicts with --input-type"):
            load_samples(path, "isolate_genomes")

    def test_old_input_type_spellings_are_rejected(self):
        for spelling in ["isolate_genome", "mag", "mags", "MAG"]:
            with self.subTest(spelling=spelling):
                path = self.manifest(f"sample_id\tinput_type\tgenome\ns\t{spelling}\tgenome.fna\n")
                with self.assertRaisesRegex(ValueError, "input_type must be one of"):
                    load_samples(path)
                with self.assertRaisesRegex(ValueError, "input_type must be one of"):
                    load_samples(path, spelling)

    def test_global_type_still_requires_unique_sample_ids(self):
        path = self.manifest("sample_id\tgenome\ns\tgenome.fna\ns\tgenome.fna\n")
        with self.assertRaisesRegex(ValueError, "duplicate sample_id"):
            load_samples(path, "MAGs")

    def test_rejects_duplicate_and_unsafe_identifiers(self):
        for records, expected in [
            ("a\tMAGs\tgenome.fna\na\tMAGs\tgenome.fna\n", "duplicate sample_id"),
            ("../escape\tMAGs\tgenome.fna\n", "Invalid sample_id"),
            ("two samples\tMAGs\tgenome.fna\n", "Invalid sample_id"),
        ]:
            with self.subTest(records=records):
                with self.assertRaisesRegex(ValueError, expected):
                    load_samples(self.manifest("sample_id\tinput_type\tgenome\n" + records))

    def test_rejects_incompatible_missing_and_empty_inputs(self):
        empty = self.root / "empty.fa"
        empty.touch()
        for sample, expected in [
            (Sample("s", "metagenome_reads", read1=self.r1), "paired"),
            (Sample("s", "isolate_reads", read1=self.r1, read2=self.r1), "different files"),
            (Sample("s", "MAGs", genome=self.genome, read1=self.r1), "cannot accompany"),
            (Sample("s", "MAGs"), "genome is required"),
            (Sample("s", "MAGs", genome=empty), "is empty"),
            (Sample("s", "MAGs", genome=self.root / "absent.fa"), "existing regular file"),
            (Sample("s", "unknown", genome=self.genome), "input_type"),
        ]:
            with self.subTest(sample=sample):
                with self.assertRaisesRegex(ValueError, expected):
                    validate_sample(sample)

    def test_rejects_malformed_manifests(self):
        for content in [
            "sample,input_type\ns,MAGs\n",
            "sample_id\tinput_type\tgenome\n",
            "sample_id\tinput_type\tgenome\tgenome\ns\tMAGs\tgenome.fna\tgenome.fna\n",
            "sample_id\tinput_type\tgenome\ns\tMAGs\tgenome.fna\textra\n",
        ]:
            with self.subTest(content=content):
                with self.assertRaises(ValueError):
                    load_samples(self.manifest(content))

    def test_missing_trailing_metadata_is_an_empty_cell(self):
        header = "sample_id\tinput_type\tgenome\tgroup\n"
        for record in ["s\tMAGs\tgenome.fna\n", "s\tMAGs\tgenome.fna\t\n"]:
            with self.subTest(record=record):
                sample = load_samples(self.manifest(header + record))[0]
                self.assertEqual(sample.metadata, {"group": ""})
        with self.assertRaisesRegex(ValueError, "genome is required"):
            load_samples(self.manifest(header + "s\tMAGs\n"))

    def test_rejects_unterminated_quoted_metadata(self):
        path = self.manifest('sample_id,input_type,genome,group\ns,MAGs,genome.fna,"case\n', "samples.csv")
        with self.assertRaisesRegex(ValueError, "malformed CSV/TSV"):
            load_samples(path)

    def test_rejects_paired_read_hardlinks_to_the_same_file(self):
        alias = self.root / "R2_alias.fastq"
        alias.hardlink_to(self.r1)
        with self.assertRaisesRegex(ValueError, "different files"):
            validate_sample(Sample("s", "metagenome_reads", read1=self.r1, read2=alias))


if __name__ == "__main__":
    unittest.main()
