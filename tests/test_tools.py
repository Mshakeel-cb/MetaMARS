"""Tool discovery and database selection without running scientific programs."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from metamars_prep.tools import discover_tools, required_tools, resolve_databases


class DiscoveryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.environment = patch.dict(os.environ, {"PATH": str(self.bin)}, clear=True)
        self.environment.start()
        self.addCleanup(self.environment.stop)

    def executable(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("#!/bin/sh\nexit 0\n")
        path.chmod(0o755)
        return str(path)

    @staticmethod
    def listing(*prefixes):
        return subprocess.CompletedProcess([], 0, json.dumps({"envs": list(map(str, prefixes))}), "")

    def test_path_tools_are_used_without_any_subprocess_probe(self):
        fastp = self.executable(self.bin / "fastp")
        spades = self.executable(self.bin / "spades")
        self.executable(self.bin / "conda")
        with patch("metamars_prep.tools.subprocess.run") as run:
            self.assertEqual(discover_tools({"fastp", "spades"}), {"fastp": [fastp], "spades": [spades]})
        run.assert_not_called()

    def test_missing_tools_use_their_own_named_environments(self):
        conda = self.executable(self.bin / "conda")
        check = self.root / "envs" / "metamars-checkm2"
        gtdb = self.root / "envs" / "metamars-gtdbtk"
        preparation = self.root / "envs" / "metamars-preparation"
        check_exe = self.executable(check / "bin" / "checkm2")
        gtdb_exe = self.executable(gtdb / "bin" / "gtdbtk")
        metaquast = self.executable(preparation / "bin" / "metaquast.py")
        fastp = self.executable(self.bin / "fastp")
        with patch("metamars_prep.tools.subprocess.run", return_value=self.listing(check, gtdb, preparation)) as run:
            commands = discover_tools({"fastp", "checkm2", "gtdbtk", "metaquast"})
        self.assertEqual(commands["fastp"], [fastp])
        for name, prefix, exe in (("checkm2", check, check_exe), ("gtdbtk", gtdb, gtdb_exe), ("metaquast", preparation, metaquast)):
            self.assertEqual(commands[name], [conda, "run", "--no-capture-output", "-p", str(prefix), exe])
        run.assert_called_once_with([conda, "env", "list", "--json"], check=True, capture_output=True, text=True, timeout=20)

    def test_mamba_is_supported_without_conda(self):
        mamba = self.executable(self.bin / "mamba")
        prefix = self.root / "metamars-gtdbtk"
        executable = self.executable(prefix / "bin" / "gtdbtk")
        with patch("metamars_prep.tools.subprocess.run", return_value=self.listing(prefix)):
            command = discover_tools({"gtdbtk"})["gtdbtk"]
        self.assertEqual(command, [mamba, "run", "-p", str(prefix), executable])

    def test_conda_exe_supports_shell_function_installations(self):
        conda = self.executable(self.root / "condabin" / "conda")
        prefix = self.root / "metamars-checkm2"
        self.executable(prefix / "bin" / "checkm2")
        with patch.dict(os.environ, {"CONDA_EXE": conda}), patch("metamars_prep.tools.subprocess.run", return_value=self.listing(prefix)):
            self.assertEqual(discover_tools({"checkm2"})["checkm2"][0], conda)

    def test_failed_conda_lookup_falls_back_to_mamba(self):
        self.executable(self.bin / "conda")
        mamba = self.executable(self.bin / "mamba")
        prefix = self.root / "metamars-checkm2"
        self.executable(prefix / "bin" / "checkm2")
        failed = subprocess.CalledProcessError(1, ["conda", "env", "list"])
        with patch("metamars_prep.tools.subprocess.run", side_effect=[failed, self.listing(prefix)]):
            self.assertEqual(discover_tools({"checkm2"})["checkm2"][0], mamba)

    def test_unrelated_env_or_nonexecutable_file_is_not_used(self):
        self.executable(self.bin / "conda")
        other = self.root / "unrelated"
        self.executable(other / "bin" / "checkm2")
        expected = self.root / "metamars-checkm2"
        (expected / "bin").mkdir(parents=True)
        (expected / "bin" / "checkm2").write_text("not executable")
        with patch("metamars_prep.tools.subprocess.run", return_value=self.listing(other, expected)):
            with self.assertRaisesRegex(ValueError, "Missing tools: checkm2"):
                discover_tools({"checkm2"})

    def test_invalid_environment_list_produces_actionable_missing_tool_error(self):
        self.executable(self.bin / "conda")
        for output in ("not JSON", "[]", '{"envs": "invalid"}'):
            with self.subTest(output=output), patch("metamars_prep.tools.subprocess.run", return_value=subprocess.CompletedProcess([], 0, output)):
                with self.assertRaisesRegex(ValueError, "metamars-checkm2"):
                    discover_tools({"checkm2"})

    def test_empty_request_and_unknown_tool(self):
        self.assertEqual(discover_tools(set()), {})
        with self.assertRaisesRegex(ValueError, "Unknown tool names: unknown"):
            discover_tools({"unknown"})

    def test_required_tools_follow_inputs_and_taxonomy_choice(self):
        genome_tools = required_tools([SimpleNamespace(input_type="isolate_genomes"), SimpleNamespace(input_type="MAGs")], True)
        self.assertEqual(genome_tools, {"checkm2", "quast"})
        reads_tools = required_tools([SimpleNamespace(input_type="metagenome_reads")])
        self.assertEqual(reads_tools, {"checkm2", "gtdbtk", "fastp", "spades", "minimap2", "samtools", "depth", "metabat2", "metaquast"})


class DatabaseTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.check = self.root / "checkm2.dmnd"
        self.check.write_bytes(b"database fixture")
        self.gtdb = self.root / "gtdb"
        self.gtdb.mkdir()
        (self.gtdb / "metadata.txt").write_text("database fixture")
        environment = patch.dict(os.environ, {}, clear=True)
        environment.start()
        self.addCleanup(environment.stop)

    def test_environment_paths_have_no_inventory_or_hash_fields(self):
        with patch.dict(os.environ, {"CHECKM2DB": str(self.check), "GTDBTK_DATA_PATH": str(self.gtdb)}):
            self.assertEqual(resolve_databases(None, None, False), {
                "checkm2": {"path": str(self.check)}, "gtdbtk": {"path": str(self.gtdb)},
            })

    def test_explicit_paths_take_precedence_and_skipped_taxonomy_needs_no_database(self):
        with patch.dict(os.environ, {"CHECKM2DB": "/missing", "GTDBTK_DATA_PATH": "/missing"}):
            self.assertEqual(resolve_databases(str(self.check), str(self.gtdb), False)["gtdbtk"]["path"], str(self.gtdb))
            self.assertEqual(resolve_databases(str(self.check), None, True), {"checkm2": {"path": str(self.check)}})

    def test_checkm_directory_requires_one_database(self):
        self.assertEqual(resolve_databases(str(self.root), None, True)["checkm2"]["path"], str(self.check))
        (self.root / "second.dmnd").write_bytes(b"other")
        with self.assertRaisesRegex(ValueError, "exactly one"):
            resolve_databases(str(self.root), None, True)

    def test_missing_or_empty_databases_fail_clearly(self):
        with self.assertRaisesRegex(ValueError, "CHECKM2DB"):
            resolve_databases(None, None, True)
        with self.assertRaisesRegex(ValueError, "GTDBTK_DATA_PATH"):
            resolve_databases(str(self.check), None, False)
        empty = self.root / "empty.dmnd"
        empty.touch()
        with self.assertRaisesRegex(ValueError, "Missing or empty CheckM2"):
            resolve_databases(str(empty), None, True)
        empty_dir = self.root / "empty"
        empty_dir.mkdir()
        with self.assertRaisesRegex(ValueError, "Missing or empty GTDB-Tk"):
            resolve_databases(str(self.check), str(empty_dir), False)


if __name__ == "__main__":
    unittest.main()
