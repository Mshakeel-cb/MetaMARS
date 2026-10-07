"""Reference downloads use isolated fixtures; these tests never access the network."""

import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

from metamars_prep.databases import ensure_checkm2
from metamars_prep.runner import StageError


class CheckM2DownloadTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.cache = self.root / "cache"
        self.runner = Mock()
        self.command = ["conda", "run", "-n", "metamars-checkm2", "checkm2"]
        environment = patch.dict(os.environ, {}, clear=True)
        environment.start()
        self.addCleanup(environment.stop)

    def database(self, path, content=b"fixture database"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        return path

    def download(self, stage, commands, cwd):
        self.assertEqual(stage, "database.checkm2")
        self.assertEqual(commands, [[
            *self.command, "database", "--download", "--path", str(cwd), "--no_write_json_db",
        ]])
        self.assertEqual(cwd.parent, self.cache)
        self.assertNotEqual(cwd, self.cache / "checkm2")
        self.database(cwd / "CheckM2_database" / "uniref.dmnd")

    def test_download_then_reuse_published_database(self):
        self.runner.run.side_effect = self.download
        first = ensure_checkm2(None, self.cache, self.runner, self.command)
        expected = self.cache / "checkm2" / "CheckM2_database" / "uniref.dmnd"
        self.assertEqual(first, {"path": str(expected)})
        self.assertEqual(expected.read_bytes(), b"fixture database")
        self.assertEqual(list(self.cache.iterdir()), [self.cache / "checkm2"])
        self.assertEqual(ensure_checkm2(None, self.cache, self.runner, self.command), first)
        self.runner.run.assert_called_once()

    def test_empty_existing_cache_can_be_filled(self):
        (self.cache / "checkm2").mkdir(parents=True)
        self.runner.run.side_effect = self.download
        result = ensure_checkm2(None, self.cache, self.runner, self.command)
        self.assertTrue(Path(result["path"]).is_file())

    def test_explicit_database_takes_precedence_over_environment(self):
        supplied = self.database(self.root / "supplied.dmnd")
        with patch.dict(os.environ, {"CHECKM2DB": "/missing"}):
            result = ensure_checkm2(str(supplied), self.cache, self.runner, self.command)
        self.assertEqual(result, {"path": str(supplied)})
        self.runner.run.assert_not_called()
        self.assertFalse(self.cache.exists())

    def test_environment_database_directory_is_reused(self):
        supplied = self.database(self.root / "supplied" / "nested" / "reference.dmnd")
        with patch.dict(os.environ, {"CHECKM2DB": str(self.root / "supplied")}):
            result = ensure_checkm2(None, self.cache, self.runner, self.command)
        self.assertEqual(result, {"path": str(supplied)})
        self.runner.run.assert_not_called()

    def test_invalid_explicit_database_does_not_download(self):
        with self.assertRaisesRegex(ValueError, "Missing or empty CheckM2"):
            ensure_checkm2(str(self.root / "missing.dmnd"), self.cache, self.runner, self.command)
        self.runner.run.assert_not_called()

    def test_invalid_environment_database_does_not_download(self):
        with patch.dict(os.environ, {"CHECKM2DB": str(self.root / "missing.dmnd")}):
            with self.assertRaisesRegex(ValueError, "Missing or empty CheckM2"):
                ensure_checkm2(None, self.cache, self.runner, self.command)
        self.runner.run.assert_not_called()

    def test_explicit_empty_and_ambiguous_databases_fail(self):
        supplied = self.root / "supplied"
        supplied.mkdir()
        with self.assertRaisesRegex(ValueError, "exactly one"):
            ensure_checkm2(str(supplied), self.cache, self.runner, self.command)
        first = self.database(supplied / "first.dmnd", b"")
        with self.assertRaisesRegex(ValueError, "Missing or empty CheckM2"):
            ensure_checkm2(str(supplied), self.cache, self.runner, self.command)
        first.write_bytes(b"first")
        self.database(supplied / "second.dmnd")
        with self.assertRaisesRegex(ValueError, "exactly one"):
            ensure_checkm2(str(supplied), self.cache, self.runner, self.command)
        self.runner.run.assert_not_called()

    def test_existing_incomplete_cache_is_preserved(self):
        partial = self.database(self.cache / "checkm2" / "partial.tar.gz")
        with self.assertRaisesRegex(ValueError, "move this directory aside"):
            ensure_checkm2(None, self.cache, self.runner, self.command)
        self.assertEqual(partial.read_bytes(), b"fixture database")
        self.runner.run.assert_not_called()

    def test_failed_download_is_not_published(self):
        def fail(stage, commands, cwd):
            self.download(stage, commands, cwd)
            raise StageError("download interrupted")
        self.runner.run.side_effect = fail
        with self.assertRaisesRegex(StageError, "download interrupted"):
            ensure_checkm2(None, self.cache, self.runner, self.command)
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_success_without_valid_database_is_not_published(self):
        for outputs in ({}, {"reference.dmnd": b""}, {"first.dmnd": b"a", "second.dmnd": b"b"}):
            with self.subTest(outputs=outputs):
                def download(stage, commands, cwd):
                    for name, data in outputs.items():
                        self.database(cwd / name, data)
                self.runner.run.side_effect = download
                with self.assertRaisesRegex(ValueError, "No database was added to the cache"):
                    ensure_checkm2(None, self.cache, self.runner, self.command)
                self.assertEqual(list(self.cache.iterdir()), [])

    def test_another_run_publishing_first_is_reused_without_overwrite(self):
        def concurrent_download(stage, commands, cwd):
            self.download(stage, commands, cwd)
            self.database(self.cache / "checkm2" / "existing.dmnd", b"other run")
        self.runner.run.side_effect = concurrent_download
        result = ensure_checkm2(None, self.cache, self.runner, self.command)
        expected = self.cache / "checkm2" / "existing.dmnd"
        self.assertEqual(result, {"path": str(expected)})
        self.assertEqual(expected.read_bytes(), b"other run")
        self.assertEqual(list(self.cache.iterdir()), [self.cache / "checkm2"])


if __name__ == "__main__":
    unittest.main()
