"""Download and reuse the official GTDB-Tk R232 reference package."""

from __future__ import annotations

import hashlib
from http.client import HTTPException
import logging
import os
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
import time
from urllib.request import Request, urlopen


# Release-specific URL and archive MD5 published in the GTDB-Tk installation guide:
# https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data
GTDB_URL = (
    "https://data.gtdb.ecogenomic.org/releases/release232/232.0/"
    "auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"
)
GTDB_MD5 = "25a59e0352b1fd150c589f56559767d4"
GTDB_DIRECTORIES = (
    "metadata", "markers", "masks", "msa", "pplacer", "radii", "taxonomy",
    "skani", "mrca_red", "split",
)
LOGGER = logging.getLogger("metamars")


def validate_gtdbtk(path: str | Path) -> dict:
    """Check the extracted R232 package root, without hashing its whole contents."""
    root = Path(path).expanduser().resolve()
    metadata = root / "metadata" / "metadata.txt"
    if not metadata.is_file():
        raise ValueError(
            f"Invalid GTDB-Tk database: {root}. Select the extracted package root "
            "containing metadata/metadata.txt, not its archive or parent directory."
        )
    try:
        settings = dict(
            line.strip().split("=", 1)
            for line in metadata.read_text(encoding="utf-8").splitlines() if "=" in line
        )
    except UnicodeError as exc:
        raise ValueError(f"Invalid GTDB-Tk metadata: {metadata}") from exc
    if settings.get("VERSION_DATA") != "r232":
        raise ValueError(f"GTDB-Tk 2.7 requires R232 reference data: {root}")
    missing = [name for name in GTDB_DIRECTORIES
               if not (root / name).is_dir() or not any((root / name).iterdir())]
    if missing:
        raise ValueError(f"Incomplete GTDB-Tk database {root}; missing/empty: {', '.join(missing)}")
    return {"path": str(root)}


def _download(destination: Path) -> None:
    LOGGER.info("Downloading GTDB-Tk R232 to %s (archive approximately 57 GB)", destination)
    request = Request(GTDB_URL, headers={"User-Agent": "metamars_prep"})
    digest = hashlib.md5(usedforsecurity=False)
    size, last_update = 0, time.monotonic()
    try:
        with urlopen(request, timeout=60) as response, destination.open("wb") as output:
            while chunk := response.read(8 * 1024 * 1024):
                output.write(chunk)
                digest.update(chunk)
                size += len(chunk)
                if time.monotonic() - last_update >= 10:
                    LOGGER.info("GTDB-Tk download: %.2f GB received", size / 1e9)
                    last_update = time.monotonic()
    except (OSError, HTTPException) as exc:
        raise ValueError(f"GTDB-Tk download failed: {exc}. Rerun to retry.") from exc
    if digest.hexdigest() != GTDB_MD5:
        raise ValueError("GTDB-Tk archive checksum mismatch; incomplete/corrupt download discarded. Rerun to retry.")
    LOGGER.info("GTDB-Tk archive downloaded (%.2f GB); checksum verified", size / 1e9)


def _safe_member(member: tarfile.TarInfo, destination: str) -> tarfile.TarInfo | None:
    if PurePosixPath(member.name).is_absolute() or ".." in PurePosixPath(member.name).parts:
        raise ValueError(f"Unsafe path in GTDB-Tk archive: {member.name}")
    # Also reject links escaping the extraction directory and special device files.
    return tarfile.data_filter(member, destination)


def ensure_gtdbtk(path: str | None, cache: Path) -> dict:
    """Use an explicit/environment path or install R232 once in the shared cache."""
    supplied = path or os.environ.get("GTDBTK_DATA_PATH")
    if supplied:
        result = validate_gtdbtk(supplied)
        LOGGER.info("Using GTDB-Tk database: %s", result["path"])
        return result
    cache = cache.expanduser().resolve()
    target = cache / "gtdbtk-r232"
    if target.exists() or target.is_symlink():
        result = validate_gtdbtk(target)
        LOGGER.info("Using cached GTDB-Tk database: %s", result["path"])
        return result
    if not hasattr(tarfile, "data_filter"):
        raise ValueError(
            "Automatic GTDB-Tk extraction requires Python with tarfile.data_filter. "
            "Update Python to a maintained security patch release, or provide --gtdbtk-db."
        )
    cache.mkdir(parents=True, exist_ok=True)
    # Only a checksum-verified, successfully extracted package gets its final name.
    with tempfile.TemporaryDirectory(prefix=".gtdbtk-r232-", dir=cache) as temporary:
        temporary = Path(temporary)
        archive = temporary / "gtdbtk_r232_data.tar.gz"
        extracted = temporary / "extracted"
        extracted.mkdir()
        _download(archive)
        LOGGER.info("Extracting GTDB-Tk R232 reference package; this may take several minutes")
        try:
            with tarfile.open(archive, mode="r:gz") as bundle:
                bundle.extractall(extracted, filter=_safe_member)
        except tarfile.TarError as exc:
            raise ValueError(f"Invalid or unsafe GTDB-Tk archive: {exc}") from exc
        # Upstream packages wrap their reference directories in a release directory.
        roots = [extracted] if (extracted / "metadata").is_dir() else [
            child for child in extracted.iterdir()
            if child.is_dir() and (child / "metadata").is_dir()
        ]
        if len(roots) != 1:
            raise ValueError("GTDB-Tk archive does not contain one extracted reference package")
        root = Path(validate_gtdbtk(roots[0])["path"])
        if target.exists():  # Another process may have completed the same download.
            return validate_gtdbtk(target)
        root.rename(target)
    LOGGER.info("GTDB-Tk R232 database ready: %s", target)
    return {"path": str(target)}
