"""Download missing reference data once and reuse a local database cache."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile


def validate_checkm2(path: str | Path) -> dict[str, str]:
    """Resolve a CheckM2 database file, including a nested download directory."""
    database = Path(path).expanduser().resolve()
    if database.is_dir():
        candidates = sorted(database.rglob("*.dmnd"))
        if len(candidates) != 1:
            raise ValueError(
                f"CheckM2 directory must contain exactly one .dmnd file: {database}. "
                "Specify the intended database file with --checkm2-db."
            )
        database = candidates[0]
    if not database.is_file() or database.stat().st_size == 0:
        raise ValueError(f"Missing or empty CheckM2 database: {database}")
    return {"path": str(database)}


def ensure_checkm2(
    path: str | None, cache: Path, runner, command: list[str],
) -> dict[str, str]:
    """Use supplied reference data or download it with CheckM2's own command.

    Supplied paths are authoritative: an invalid path is an error, not a reason
    to download another copy. Downloads are validated in a temporary directory
    before publication, without changing CheckM2's installed configuration.
    """
    supplied = path or os.environ.get("CHECKM2DB")
    if supplied:
        return validate_checkm2(supplied)

    cache = Path(cache).expanduser().resolve()
    destination = cache / "checkm2"
    if destination.exists():
        try:
            return validate_checkm2(destination)
        except ValueError as error:
            if not destination.is_dir() or any(destination.iterdir()):
                raise ValueError(
                    f"Existing CheckM2 cache is incomplete or ambiguous: {destination}. "
                    "Supply a valid --checkm2-db path, or move this directory aside "
                    "before downloading again. Existing files were not changed."
                ) from error

    cache.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".checkm2-download-", dir=cache) as directory:
        temporary = Path(directory)
        runner.run(
            "database.checkm2",
            [[*command, "database", "--download", "--path", str(temporary), "--no_write_json_db"]],
            temporary,
        )
        try:
            downloaded = Path(validate_checkm2(temporary)["path"])
        except ValueError as error:
            raise ValueError(
                "CheckM2 download did not produce exactly one nonempty .dmnd file. "
                "See logs/database.checkm2.log. No database was added to the cache."
            ) from error
        relative = downloaded.relative_to(temporary)
        # Remove only an empty destination; never overwrite an existing download.
        if destination.exists():
            if destination.is_dir() and not any(destination.iterdir()):
                destination.rmdir()
            else:
                # Another run may have finished downloading while this one ran.
                return validate_checkm2(destination)
        temporary.rename(destination)
    return {"path": str(destination / relative)}
