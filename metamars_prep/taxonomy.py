"""Read genome quality and taxonomy reports and create lineage comparison groups.

Genome identifiers are opaque keys: never trim, case-fold, or shorten them.
Quality exclusions and unresolved classifications remain in the returned tables.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterator


RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
_PREFIXES = dict(zip(("d", "p", "c", "o", "f", "g", "s"), RANKS))


def _report_rows(path: Path, required: set[str]) -> Iterator[tuple[int, dict[str, str]]]:
    """Read a TSV with a validated header and structurally complete rows."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        try:
            columns = reader.fieldnames
            if not columns:
                raise ValueError(f"{path}: missing report header")
            if len(columns) != len(set(columns)):
                raise ValueError(f"{path}: duplicate report columns")
            missing = required.difference(columns)
            if missing:
                raise ValueError(f"{path}: missing required columns: {', '.join(sorted(missing))}")
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"{path}:{reader.line_num}: row does not match report columns")
                yield reader.line_num, row
        except csv.Error as error:
            raise ValueError(f"{path}:{reader.line_num}: malformed TSV: {error}") from error


def _identifier(value: str, path: Path, line: int, column: str) -> str:
    if not value.strip():
        raise ValueError(f"{path}:{line}: empty {column}")
    return value


def _finite_number(value: str, path: Path, line: int, column: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{path}:{line}: invalid {column}: {value!r}") from error
    if not math.isfinite(number):
        raise ValueError(f"{path}:{line}: nonfinite {column}: {value!r}")
    return number


def parse_checkm2(path: Path) -> dict[str, dict]:
    """Index a CheckM2 quality report by exact Name, retaining the source row."""
    results: dict[str, dict] = {}
    for line, row in _report_rows(path, {"Name", "Completeness", "Contamination"}):
        genome_id = _identifier(row["Name"], path, line, "Name")
        if genome_id in results:
            raise ValueError(f"{path}:{line}: duplicate CheckM2 Name: {genome_id!r}")
        completeness = _finite_number(row["Completeness"], path, line, "Completeness")
        contamination = _finite_number(row["Contamination"], path, line, "Contamination")
        if not 0 <= completeness <= 100:
            raise ValueError(f"{path}:{line}: Completeness must be between 0 and 100: {completeness}")
        if contamination < 0:
            raise ValueError(f"{path}:{line}: Contamination must be nonnegative: {contamination}")
        results[genome_id] = {
            "completeness": completeness,
            "contamination": contamination,
            "raw": row,
        }
    return results


def _lineage(classification: str) -> dict[str, str]:
    ranks = dict.fromkeys(RANKS, "")
    for component in classification.split(";"):
        component = component.strip()
        prefix, separator, name = component.partition("__")
        if not separator or prefix not in _PREFIXES:
            continue
        rank = _PREFIXES[prefix]
        if not name.strip():
            continue
        if ranks[rank] and ranks[rank] != component:
            raise ValueError(f"conflicting {rank} assignments in classification: {classification!r}")
        ranks[rank] = component
    return ranks


def parse_gtdb(paths: list[Path]) -> dict[str, dict]:
    """Combine GTDB-Tk summaries, rejecting conflicting duplicate genome rows.

    Identical repeated rows are accepted. Unclassified ranks are empty strings;
    the original classification, warnings and other report columns are retained.
    """
    results: dict[str, dict] = {}
    for path in paths:
        for line, row in _report_rows(path, {"user_genome", "classification"}):
            genome_id = _identifier(row["user_genome"], path, line, "user_genome")
            try:
                ranks = _lineage(row["classification"])
            except ValueError as error:
                raise ValueError(f"{path}:{line}: {error}") from error
            record = {
                **row,
                **ranks,
                "warnings": row.get("warnings", ""),
                "classification_method": row.get("classification_method", ""),
                "raw": row,
            }
            if genome_id in results and results[genome_id] != record:
                raise ValueError(f"{path}:{line}: conflicting GTDB-Tk genome: {genome_id!r}")
            results[genome_id] = record
    return results


def collect_gtdb_failures(report_dir: Path) -> dict[str, str]:
    """Collect documented, headerless GTDB-Tk genome exclusion reports.

    Only reports produced with the default ``gtdbtk`` prefix are recognized.
    Each row is the exact genome identifier and reason separated by one tab.
    Root-level symlinks and repeated reasons are deduplicated. Distinct reasons
    for the same genome are retained in deterministic order, joined with ``; ``.
    """
    names = {
        "gtdbtk.failed_genomes.tsv",
        "gtdbtk.bac120.filtered.tsv",
        "gtdbtk.ar53.filtered.tsv",
        "gtdbtk.align.failed.tsv",
    }
    report_dir = Path(report_dir)
    if not report_dir.is_dir():
        raise ValueError(f"Missing GTDB-Tk report directory: {report_dir}")
    seen: set[Path] = set()
    reasons: dict[str, set[str]] = {}
    for path in sorted(report_dir.rglob("*.tsv")):
        if path.name not in names:
            continue
        resolved = path.resolve(strict=True)
        if resolved in seen:
            continue
        seen.add(resolved)
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for line_number, line in enumerate(handle, 1):
                text = line.rstrip("\r\n")
                if not text:
                    continue
                fields = text.split("\t")
                if len(fields) != 2:
                    raise ValueError(f"{path}:{line_number}: expected genome ID and failure reason separated by a tab")
                genome_id = _identifier(fields[0], path, line_number, "genome ID")
                reason = fields[1]
                if not reason.strip():
                    raise ValueError(f"{path}:{line_number}: empty GTDB-Tk failure reason")
                reasons.setdefault(genome_id, set()).add(reason)
    return {genome_id: "; ".join(sorted(values)) for genome_id, values in sorted(reasons.items())}


def _eligible(value: object) -> bool:
    if isinstance(value, str):
        normalized = value.strip().casefold()
        if normalized in {"true", "yes", "y", "1", "on"}:
            return True
        if normalized in {"false", "no", "n", "0", "off", ""}:
            return False
        raise ValueError(f"invalid comparison_eligible boolean: {value!r}")
    return bool(value)


def group_genomes(genomes: list[dict], taxonomy: dict[str, dict], rank: str) -> list[dict]:
    """Return a comparison membership row for every input genome.

    Only eligible genomes assigned at the selected rank receive a group ID.
    IDs depend on all ancestral rank assignments, keeping homonyms separate.
    """
    if rank not in RANKS:
        raise ValueError(f"invalid taxonomic rank {rank!r}; choose from {', '.join(RANKS)}")
    selected_ranks = RANKS[: RANKS.index(rank) + 1]
    groups: list[dict] = []
    for genome in genomes:
        genome_id = genome["genome_id"]
        tax = taxonomy.get(genome_id, {})
        lineage = _lineage(tax.get("classification", ""))
        for tax_rank in RANKS:
            value = tax.get(tax_rank, lineage[tax_rank])
            if value and "__" in value and value.partition("__")[2].strip():
                lineage[tax_rank] = value
        taxon = lineage[rank]
        group_id = ""
        if not _eligible(genome.get("comparison_eligible", False)):
            status = "quality_excluded"
        elif not taxon:
            status = "unresolved"
        else:
            status = "assigned"
            key = json.dumps([lineage[name] for name in selected_ranks], ensure_ascii=False, separators=(",", ":"))
            digest = hashlib.sha256(key.encode("utf-8")).hexdigest()[:20]
            group_id = f"{rank}_{digest}"
        groups.append({
            "sample_id": genome["sample_id"],
            "genome_id": genome_id,
            "tax_rank": rank,
            "taxon": taxon,
            "group_id": group_id,
            "status": status,
        })
    return groups
