import gzip
import hashlib
import tempfile
import unittest
from pathlib import Path

from metamars_prep.registry import build_membership, iter_fasta, normalize_fasta, write_tsv


class RegistryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def fasta(self, name, contents):
        path = self.root / name
        path.write_text(contents)
        return path

    def test_normalization_preserves_sequence_identity_and_distinguishes_samples(self):
        source = self.root / "original.fasta.gz"
        with gzip.open(source, "wt") as handle:
            handle.write(">old description\nacgtnry\n>second\nGGCC\n")
        dest = self.root / "canonical.fa.gz"
        rows = normalize_fasta(source, dest, "sample1")
        self.assertEqual(list(iter_fasta(dest)), [
            ("sample1__asm1__contig000001", "ACGTNRY"),
            ("sample1__asm1__contig000002", "GGCC"),
        ])
        self.assertEqual(rows[0]["original_id"], "old")
        self.assertEqual(rows[0]["sequence_sha256"], hashlib.sha256(b"ACGTNRY").hexdigest())
        self.assertEqual(rows[0]["length"], 7)
        self.assertEqual(rows[1]["gc_percent"], 100)
        second = normalize_fasta(source, self.root / "other.fa", "sample2")
        self.assertNotEqual(rows[0]["contig_id"], second[0]["contig_id"])
        self.assertEqual(rows[0]["sequence_sha256"], second[0]["sequence_sha256"])

    def test_membership_retains_short_unbinned_contigs_without_synthetic_hosts(self):
        assembly = self.fasta("assembly.fa", ">a\nACGTACGT\n>b\nNN\n>c\nGGCCAATT\n")
        bin_path = self.fasta("bin.fa", ">a\nACGTACGT\n")
        output = self.root / "unbinned.fa"
        rows = build_membership(assembly, {"sample_bin1": bin_path}, output, 5)
        self.assertEqual(list(iter_fasta(output)), [("b", "NN"), ("c", "GGCCAATT")])
        self.assertEqual([row["binning_status"] for row in rows], ["binned", "too_short", "unbinned"])
        self.assertEqual([row["genome_id"] for row in rows], ["sample_bin1", "", ""])
        self.assertEqual([row["host_assignment_status"] for row in rows], ["bin_supported", "unresolved", "unresolved"])
        self.assertEqual(len(rows), 3)

    def test_identical_sequences_keep_distinct_contig_and_bin_identities(self):
        assembly = self.fasta("assembly.fa", ">copy1\nACGT\n>copy2\nACGT\n")
        bin_path = self.fasta("bin.fa", ">copy1\nACGT\n")
        output = self.root / "unbinned.fa"
        rows = build_membership(assembly, {"genome1": bin_path}, output, 1)
        self.assertEqual([(row["contig_id"], row["genome_id"]) for row in rows], [("copy1", "genome1"), ("copy2", "")])
        self.assertEqual(list(iter_fasta(output)), [("copy2", "ACGT")])

    def test_all_binned_and_no_bins_are_valid_boundaries(self):
        assembly = self.fasta("assembly.fa", ">a\nACGT\n")
        output = self.root / "unbinned.fa"
        rows = build_membership(assembly, {"g": assembly}, output, 5)
        self.assertEqual(rows[0]["binning_status"], "binned")
        self.assertEqual(output.read_text(), "")
        rows = build_membership(assembly, {}, output, 5)
        self.assertEqual(rows[0]["binning_status"], "too_short")
        self.assertEqual(list(iter_fasta(output)), [("a", "ACGT")])

    def test_rejects_foreign_changed_and_multiply_assigned_contigs(self):
        assembly = self.fasta("assembly.fa", ">a\nACGT\n")
        valid = self.fasta("valid.fa", ">a\nACGT\n")
        foreign = self.fasta("foreign.fa", ">z\nACGT\n")
        changed = self.fasta("changed.fa", ">a\nTGCA\n")
        output = self.fasta("unbinned.fa", "previous result\n")
        for bins, expected in [
            ({"g": foreign}, "absent"),
            ({"g": changed}, "differs"),
            ({"g": valid, "h": valid}, "multiple bins"),
        ]:
            with self.subTest(bins=bins):
                with self.assertRaisesRegex(ValueError, expected):
                    build_membership(assembly, bins, output, 5)
                self.assertEqual(output.read_text(), "previous result\n")

    def test_rejects_invalid_fasta_and_preserves_previous_output(self):
        dest = self.fasta("canonical.fa", ">previous\nACGT\n")
        for contents in ["", "ACGT\n", ">a\nACGU\n", ">a\n", ">\nACGT\n", ">a\nACGT\n>a\nTTTT\n"]:
            with self.subTest(contents=contents):
                source = self.fasta("invalid.fa", contents)
                with self.assertRaises(ValueError):
                    normalize_fasta(source, dest, "sample1")
                self.assertEqual(dest.read_text(), ">previous\nACGT\n")

    def test_write_tsv_has_stable_header_for_empty_results(self):
        target = self.root / "empty.tsv"
        write_tsv(target, [], ["contig_id", "genome_id"])
        self.assertEqual(target.read_text(), "contig_id\tgenome_id\n")
        write_tsv(target, [{"contig_id": "a", "genome_id": "g", "raw_metadata": "retained elsewhere"}], ["contig_id", "genome_id"])
        self.assertEqual(target.read_text(), "contig_id\tgenome_id\na\tg\n")


if __name__ == "__main__":
    unittest.main()
