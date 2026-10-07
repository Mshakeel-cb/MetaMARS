"""Offline tests for database selection, verified downloads and safe extraction."""

import hashlib
import io
import os
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch
from urllib.error import URLError

from metamars_prep.gtdb_download import GTDB_DIRECTORIES, ensure_gtdbtk, validate_gtdbtk


class GTDBDownloadTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.cache = self.root / "cache"
        environment = patch.dict(os.environ, {}, clear=True)
        environment.start()
        self.addCleanup(environment.stop)

    def reference(self, directory):
        for name in GTDB_DIRECTORIES:
            folder = directory / name
            folder.mkdir(parents=True)
            (folder / "data").write_text("reference")
        (directory / "metadata" / "metadata.txt").write_text("VERSION_DATA=r232\n")
        return directory

    def archive(self, extra=None, version="r232", wrapped=True):
        data = io.BytesIO()
        prefix = "release232/" if wrapped else ""
        with tarfile.open(fileobj=data, mode="w:gz") as archive:
            for name in GTDB_DIRECTORIES:
                path = prefix + name + ("/metadata.txt" if name == "metadata" else "/data")
                content = f"VERSION_DATA={version}\n".encode() if name == "metadata" else b"reference"
                member = tarfile.TarInfo(path)
                member.size = len(content)
                archive.addfile(member, io.BytesIO(content))
            if extra is not None:
                archive.addfile(extra, io.BytesIO(b"x") if extra.isfile() else None)
        return data.getvalue()

    def download(self, data):
        md5 = hashlib.md5(data, usedforsecurity=False).hexdigest()
        checksum = patch("metamars_prep.gtdb_download.GTDB_MD5", md5)
        checksum.start()
        self.addCleanup(checksum.stop)
        request = patch("metamars_prep.gtdb_download.urlopen", return_value=io.BytesIO(data))
        opened = request.start()
        self.addCleanup(request.stop)
        return opened

    def test_explicit_database_does_not_download_or_create_cache(self):
        reference = self.reference(self.root / "existing")
        with patch("metamars_prep.gtdb_download.urlopen") as opened:
            self.assertEqual(ensure_gtdbtk(str(reference), self.cache), {"path": str(reference)})
        opened.assert_not_called()
        self.assertFalse(self.cache.exists())

    def test_environment_database_and_cli_precedence(self):
        reference = self.reference(self.root / "existing")
        with patch.dict(os.environ, {"GTDBTK_DATA_PATH": str(reference)}):
            self.assertEqual(ensure_gtdbtk(None, self.cache)["path"], str(reference))
            with self.assertRaisesRegex(ValueError, "Invalid GTDB-Tk"):
                ensure_gtdbtk(str(self.root / "missing"), self.cache)

    def test_invalid_supplied_database_does_not_trigger_download(self):
        with patch("metamars_prep.gtdb_download.urlopen") as opened:
            with self.assertRaisesRegex(ValueError, "extracted package root"):
                ensure_gtdbtk(str(self.root), self.cache)
        opened.assert_not_called()

    def test_validate_rejects_empty_component_and_wrong_release(self):
        reference = self.reference(self.root / "existing")
        (reference / "skani" / "data").unlink()
        with self.assertRaisesRegex(ValueError, "skani"):
            validate_gtdbtk(reference)
        (reference / "metadata" / "metadata.txt").write_text("VERSION_DATA=r214\n")
        with self.assertRaisesRegex(ValueError, "requires R232"):
            validate_gtdbtk(reference)

    def test_cache_is_reused_without_network(self):
        reference = self.reference(self.cache / "gtdbtk-r232")
        with patch("metamars_prep.gtdb_download.urlopen") as opened:
            self.assertEqual(ensure_gtdbtk(None, self.cache)["path"], str(reference))
        opened.assert_not_called()

    def test_download_verifies_extracts_and_reuses_wrapped_package(self):
        opened = self.download(self.archive())
        with self.assertLogs("metamars", level="INFO") as logs:
            result = ensure_gtdbtk(None, self.cache)
        reference = self.cache / "gtdbtk-r232"
        self.assertEqual(result, {"path": str(reference)})
        self.assertEqual(validate_gtdbtk(reference), result)
        self.assertEqual(list(self.cache.iterdir()), [reference])
        self.assertIn("checksum verified", " ".join(logs.output))
        self.assertEqual(ensure_gtdbtk(None, self.cache), result)
        opened.assert_called_once()
        self.assertEqual(opened.call_args.kwargs["timeout"], 60)

    def test_unwrapped_archive_is_supported(self):
        self.download(self.archive(wrapped=False))
        self.assertTrue(Path(ensure_gtdbtk(None, self.cache)["path"]).is_dir())

    def test_bad_checksum_is_not_published(self):
        self.download(self.archive())
        with patch("metamars_prep.gtdb_download.GTDB_MD5", "bad-checksum"):
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                ensure_gtdbtk(None, self.cache)
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_download_failure_is_not_published(self):
        with patch("metamars_prep.gtdb_download.urlopen", side_effect=URLError("offline")):
            with self.assertRaisesRegex(ValueError, "download failed"):
                ensure_gtdbtk(None, self.cache)
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_partial_response_failure_is_cleaned(self):
        response = io.BytesIO()
        with patch("metamars_prep.gtdb_download.urlopen", return_value=response):
            with patch.object(response, "read", side_effect=[b"partial", OSError("lost connection")]):
                with self.assertRaisesRegex(ValueError, "download failed"):
                    ensure_gtdbtk(None, self.cache)
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_parent_traversal_is_rejected(self):
        member = tarfile.TarInfo("../../escaped")
        member.size = 1
        self.download(self.archive(extra=member))
        with self.assertRaisesRegex(ValueError, "Unsafe path"):
            ensure_gtdbtk(None, self.cache)
        self.assertFalse((self.cache / "escaped").exists())
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_escaping_symbolic_link_is_rejected(self):
        member = tarfile.TarInfo("release232/outside")
        member.type = tarfile.SYMTYPE
        member.linkname = "../../../../outside"
        self.download(self.archive(extra=member))
        with self.assertRaisesRegex(ValueError, "unsafe GTDB-Tk archive"):
            ensure_gtdbtk(None, self.cache)
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_wrong_archive_release_is_not_published(self):
        self.download(self.archive(version="r214"))
        with self.assertRaisesRegex(ValueError, "requires R232"):
            ensure_gtdbtk(None, self.cache)
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_missing_safe_extraction_support_fails_before_download(self):
        with patch("metamars_prep.gtdb_download.tarfile") as tar, patch("metamars_prep.gtdb_download.urlopen") as opened:
            del tar.data_filter
            with self.assertRaisesRegex(ValueError, "security patch"):
                ensure_gtdbtk(None, self.cache)
        opened.assert_not_called()


if __name__ == "__main__":
    unittest.main()
