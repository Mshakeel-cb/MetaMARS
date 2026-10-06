"""Report validation and lineage grouping behavior."""

import tempfile
import unittest
from pathlib import Path

from metamars_prep.taxonomy import RANKS, collect_gtdb_failures, group_genomes, parse_checkm2, parse_gtdb


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.root = Path(self.tempdir.name)

    def report(self, text, name="report.tsv"):
        path = self.root / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_checkm2_preserves_exact_ids_and_raw_fields(self):
        path = self.report("Name\tCompleteness\tContamination\tModel\n g.fa \t91.2\t2\tgeneral\ng.fa\t80\t0\tspecific\n")
        records = parse_checkm2(path)
        self.assertEqual(set(records), {" g.fa ", "g.fa"})
        self.assertEqual(records[" g.fa "]["completeness"], 91.2)
        self.assertEqual(records[" g.fa "]["contamination"], 2.0)
        self.assertEqual(records[" g.fa "]["raw"]["Model"], "general")

    def test_checkm2_rejects_duplicate_names(self):
        path = self.report("Name\tCompleteness\tContamination\ng\t90\t1\ng\t90\t1\n")
        with self.assertRaisesRegex(ValueError, "duplicate CheckM2"):
            parse_checkm2(path)

    def test_checkm2_rejects_missing_columns(self):
        with self.assertRaisesRegex(ValueError, "missing required columns: Contamination"):
            parse_checkm2(self.report("Name\tCompleteness\ng\t90\n"))

    def test_checkm2_rejects_invalid_and_nonfinite_numbers(self):
        for field in ("Completeness", "Contamination"):
            for bad in ("NaN", "inf", "-Infinity", "", "missing"):
                with self.subTest(field=field, bad=bad):
                    values = {"Completeness": "90", "Contamination": "1", field: bad}
                    path = self.report(f"Name\tCompleteness\tContamination\ng\t{values['Completeness']}\t{values['Contamination']}\n")
                    with self.assertRaisesRegex(ValueError, field):
                        parse_checkm2(path)

    def test_checkm2_rejects_invalid_ranges_but_allows_contamination_over_100(self):
        for completeness, contamination in (("-1", "0"), ("100.1", "0"), ("95", "-0.1")):
            with self.subTest(completeness=completeness, contamination=contamination):
                path = self.report(f"Name\tCompleteness\tContamination\ng\t{completeness}\t{contamination}\n")
                with self.assertRaises(ValueError):
                    parse_checkm2(path)
        path = self.report("Name\tCompleteness\tContamination\ng\t100\t120.5\n")
        self.assertEqual(parse_checkm2(path)["g"]["contamination"], 120.5)

    def test_reports_reject_duplicate_headers_and_malformed_rows(self):
        for text in (
            "Name\tCompleteness\tContamination\tName\ng\t90\t1\tg\n",
            "Name\tCompleteness\tContamination\ng\t90\n",
            "Name\tCompleteness\tContamination\ng\t90\t1\textra\n",
            "Name\tCompleteness\tContamination\n \t90\t1\n",
        ):
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_checkm2(self.report(text))

    def test_gtdb_combines_domains_and_preserves_ids_metadata_and_ranks(self):
        header = "user_genome\tclassification\twarnings\tclassification_method\n"
        bac = self.report(header + " g.fa \td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G one\tlow coverage\tANI\n", "bac.tsv")
        arc = self.report(header + "g.fa\td__Archaea;p__A;c__;o__;f__;g__;s__\t\tplacement\n", "arc.tsv")
        records = parse_gtdb([bac, arc])
        self.assertEqual(set(records), {" g.fa ", "g.fa"})
        self.assertTrue(set(RANKS).issubset(records[" g.fa "]))
        self.assertEqual(records[" g.fa "]["species"], "s__G one")
        self.assertEqual(records[" g.fa "]["warnings"], "low coverage")
        self.assertEqual(records[" g.fa "]["classification_method"], "ANI")
        self.assertEqual(records["g.fa"]["genus"], "")

    def test_gtdb_allows_identical_but_rejects_conflicting_duplicates(self):
        first = self.report("user_genome\tclassification\ng\td__Bacteria;g__G\n", "a.tsv")
        same = self.report(first.read_text(), "b.tsv")
        other = self.report("user_genome\tclassification\ng\td__Bacteria;g__H\n", "c.tsv")
        self.assertEqual(len(parse_gtdb([first, same])), 1)
        with self.assertRaisesRegex(ValueError, "conflicting GTDB-Tk"):
            parse_gtdb([first, other])

    def test_gtdb_requires_identity_and_classification_columns(self):
        with self.assertRaisesRegex(ValueError, "missing required columns: classification"):
            parse_gtdb([self.report("user_genome\ng\n")])

    def test_gtdb_rejects_conflicting_rank_components(self):
        with self.assertRaisesRegex(ValueError, "conflicting genus"):
            parse_gtdb([self.report("user_genome\tclassification\ng\td__Bacteria;g__G;g__H\n")])

    def test_gtdb_rejects_truncated_quoted_classification(self):
        path = self.report('user_genome\tclassification\ng\t"d__Bacteria;g__G\n')
        with self.assertRaisesRegex(ValueError, "malformed TSV"):
            parse_gtdb([path])

    def test_gtdb_explicit_unclassified_record_remains_a_valid_summary(self):
        path = self.report("user_genome\tclassification\tgtdb_warning\ng\tUnclassified Bacteria\tno markers\n")
        record = parse_gtdb([path])["g"]
        self.assertEqual(record["classification"], "Unclassified Bacteria")
        self.assertEqual(record["domain"], "")

    def test_failure_reports_preserve_ids_and_reasons_and_deduplicate_symlinks(self):
        identify = self.root / "identify"
        align = self.root / "align"
        identify.mkdir()
        (align / "intermediate_results").mkdir(parents=True)
        failed = identify / "gtdbtk.failed_genomes.tsv"
        failed.write_text(" g.fa \tNo genes were called by Prodigal\ng.fa\tEmpty file\n")
        (self.root / failed.name).symlink_to(failed)
        filtered = "Insufficient number of amino acids in MSA (1.2%)"
        (align / "gtdbtk.bac120.filtered.tsv").write_text(f"g2\t{filtered}\ng2\t{filtered}\n")
        (align / "gtdbtk.ar53.filtered.tsv").write_text("")
        (align / "intermediate_results" / "gtdbtk.align.failed.tsv").write_text("g3\tNo bacterial or archaeal marker\n")
        failures = collect_gtdb_failures(self.root)
        self.assertEqual(failures, {
            " g.fa ": "No genes were called by Prodigal",
            "g.fa": "Empty file",
            "g2": filtered,
            "g3": "No bacterial or archaeal marker",
        })

    def test_failure_reports_merge_distinct_reasons_and_ignore_unrecognized_names(self):
        self.report("g\treason B\ng\treason A\n", "gtdbtk.failed_genomes.tsv")
        self.report("g\treason B\n", "gtdbtk.bac120.filtered.tsv")
        for name in ("other.failed_genomes.tsv", "gtdbtk.unrelated.filtered.tsv", "gtdbtk.bac120.summary.tsv"):
            self.report("not a valid failure row", name)
        self.assertEqual(collect_gtdb_failures(self.root), {"g": "reason A; reason B"})

    def test_failure_reports_reject_malformed_rows_and_blank_identifiers_or_reasons(self):
        for text in ("g\n", "g\treason\textra\n", " \treason\n", "g\t \n"):
            with self.subTest(text=text):
                self.report(text, "gtdbtk.failed_genomes.tsv")
                with self.assertRaises(ValueError):
                    collect_gtdb_failures(self.root)

    def test_absence_of_failure_files_is_empty_and_missing_directory_fails(self):
        self.assertEqual(collect_gtdb_failures(self.root), {})
        with self.assertRaisesRegex(ValueError, "Missing GTDB-Tk"):
            collect_gtdb_failures(self.root / "missing")


class GroupingTests(unittest.TestCase):
    @staticmethod
    def genome(genome_id, eligible=True, sample_id="sample"):
        return {"sample_id": sample_id, "genome_id": genome_id, "comparison_eligible": eligible}

    @staticmethod
    def taxonomy(species="one", family="F"):
        return {"classification": f"d__Bacteria;p__P;c__C;o__O;f__{family};g__G;s__G {species}"}

    def test_genus_groups_species_while_species_keeps_them_separate(self):
        genomes = [self.genome("a", sample_id="s1"), self.genome("b", " TRUE ", "s2")]
        taxonomy = {"a": self.taxonomy("one"), "b": self.taxonomy("two")}
        genus = group_genomes(genomes, taxonomy, "genus")
        species = group_genomes(genomes, taxonomy, "species")
        self.assertEqual(genus[0]["group_id"], genus[1]["group_id"])
        self.assertNotEqual(species[0]["group_id"], species[1]["group_id"])
        self.assertEqual(genus[0]["taxon"], "g__G")
        self.assertEqual(genus[1]["sample_id"], "s2")
        self.assertEqual(genus[1]["status"], "assigned")

    def test_homonyms_have_distinct_stable_group_ids(self):
        genomes = [self.genome("a"), self.genome("b")]
        taxonomy = {"a": self.taxonomy(family="F"), "b": self.taxonomy(family="Other")}
        forward = group_genomes(genomes, taxonomy, "genus")
        reverse = group_genomes(list(reversed(genomes)), taxonomy, "genus")
        self.assertNotEqual(forward[0]["group_id"], forward[1]["group_id"])
        self.assertEqual(forward[0]["group_id"], reverse[1]["group_id"])

    def test_unresolved_rows_never_receive_group_ids(self):
        genomes = [self.genome("a"), self.genome("b"), self.genome("c")]
        taxonomy = {"a": {"classification": "d__Bacteria;g__G;s__"}, "c": {"classification": "Unclassified Bacteria"}}
        rows = group_genomes(genomes, taxonomy, "species")
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(row["status"] == "unresolved" and row["group_id"] == "" for row in rows))

    def test_quality_exclusions_remain_even_when_classified(self):
        values = [False, "False", " false ", "0", "no", "off", None, ""]
        genomes = [self.genome(str(index), value) for index, value in enumerate(values)]
        taxonomy = {genome["genome_id"]: self.taxonomy() for genome in genomes}
        rows = group_genomes(genomes, taxonomy, "species")
        self.assertEqual(len(rows), len(values))
        self.assertTrue(all(row["status"] == "quality_excluded" and row["group_id"] == "" for row in rows))
        self.assertTrue(all(row["taxon"] == "s__G one" for row in rows))

    def test_identity_lookup_does_not_trim_or_remove_extensions(self):
        genomes = [self.genome(" g.fa "), self.genome("g.fa"), self.genome("g")]
        rows = group_genomes(genomes, {" g.fa ": self.taxonomy()}, "species")
        self.assertEqual([row["status"] for row in rows], ["assigned", "unresolved", "unresolved"])
        self.assertEqual(rows[0]["genome_id"], " g.fa ")

    def test_bad_rank_and_ambiguous_boolean_fail(self):
        with self.assertRaisesRegex(ValueError, "invalid taxonomic rank"):
            group_genomes([], {}, "strain")
        with self.assertRaisesRegex(ValueError, "comparison_eligible"):
            group_genomes([self.genome("g", "maybe")], {}, "species")


if __name__ == "__main__":
    unittest.main()
