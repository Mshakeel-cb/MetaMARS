"""Canonical contig identities and explicit, lossless bin membership."""

from __future__ import annotations

import csv
import gzip
import hashlib
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

from .models import validate_identifier


DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
CONTIG_FIELDS = ["sample_id", "assembly_id", "contig_id", "original_id", "length", "gc_percent", "sequence_sha256"]
MEMBERSHIP_FIELDS = ["contig_id", "genome_id", "binning_status", "host_assignment_status"]


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Stream nonempty FASTA records as uppercase IUPAC DNA with unique IDs."""
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    seen: set[str] = set()
    identifier: str | None = None
    sequence: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    if not sequence:
                        raise ValueError(f"{path}: empty FASTA record {identifier!r}")
                    yield identifier, "".join(sequence)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA identifier")
                identifier = header.split()[0]
                if identifier in seen:
                    raise ValueError(f"{path}:{line_number}: duplicate FASTA identifier {identifier!r}")
                seen.add(identifier)
                sequence = []
            else:
                if identifier is None:
                    raise ValueError(f"{path}:{line_number}: sequence before first FASTA header")
                bases = line.upper()
                invalid = set(bases) - DNA_ALPHABET
                if invalid:
                    raise ValueError(f"{path}:{line_number}: invalid IUPAC DNA characters {''.join(sorted(invalid))!r}")
                sequence.append(bases)
    if identifier is None:
        raise ValueError(f"{path}: FASTA contains no records")
    if not sequence:
        raise ValueError(f"{path}: empty FASTA record {identifier!r}")
    yield identifier, "".join(sequence)


@contextmanager
def _atomic_text_output(path: Path) -> Iterator[TextIO]:
    """Publish complete outputs only, including when input validation fails."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        with opener(temporary, "wt", encoding="utf-8", newline="") as handle:
            yield handle
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _write_record(handle: TextIO, identifier: str, sequence: str) -> None:
    handle.write(f">{identifier}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start:start + 80] + "\n")


def _sequence_digest(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def normalize_fasta(source: Path, dest: Path, sample_id: str) -> list[dict]:
    """Create one canonical assembly; GC percentage includes ambiguous bases."""
    validate_identifier(sample_id)
    if Path(source).resolve() == Path(dest).resolve():
        raise ValueError("Source and canonical FASTA destination must differ")
    assembly_id = f"{sample_id}__asm1"
    rows = []
    with _atomic_text_output(dest) as handle:
        for index, (original_id, sequence) in enumerate(iter_fasta(source), 1):
            contig_id = f"{assembly_id}__contig{index:06d}"
            _write_record(handle, contig_id, sequence)
            rows.append({
                "sample_id": sample_id,
                "assembly_id": assembly_id,
                "contig_id": contig_id,
                "original_id": original_id,
                "length": len(sequence),
                "gc_percent": round(100 * (sequence.count("G") + sequence.count("C")) / len(sequence), 6),
                "sequence_sha256": _sequence_digest(sequence),
            })
    return rows


def build_membership(
    contigs: Path,
    bin_fastas: dict[str, Path],
    unbinned_output: Path,
    min_bin_length: int,
) -> list[dict]:
    """Retain every contig, treating bins as membership views of the assembly.

    Sequence length and SHA-256 validate bin identity without storing entire
    assemblies in memory. Short unbinned sequences remain in the output FASTA.
    """
    if min_bin_length < 0:
        raise ValueError("min_bin_length must be nonnegative")
    output_path = Path(unbinned_output).resolve()
    if output_path in {Path(contigs).resolve(), *(Path(path).resolve() for path in bin_fastas.values())}:
        raise ValueError("Unbinned output must differ from assembly and bin inputs")
    assignments: dict[str, tuple[str, int, str]] = {}
    for genome_id, fasta in bin_fastas.items():
        validate_identifier(genome_id, "genome_id")
        for contig_id, sequence in iter_fasta(fasta):
            if contig_id in assignments:
                raise ValueError(f"Contig {contig_id!r} belongs to multiple bins: {assignments[contig_id][0]!r}, {genome_id!r}")
            assignments[contig_id] = (genome_id, len(sequence), _sequence_digest(sequence))
    rows = []
    with _atomic_text_output(unbinned_output) as handle:
        for contig_id, sequence in iter_fasta(contigs):
            membership = assignments.pop(contig_id, None)
            if membership is not None:
                genome_id, length, digest = membership
                if len(sequence) != length or _sequence_digest(sequence) != digest:
                    raise ValueError(f"Bin {genome_id!r}: sequence for {contig_id!r} differs from canonical assembly")
                status, host_status = "binned", "bin_supported"
            else:
                genome_id = ""
                status = "too_short" if len(sequence) < min_bin_length else "unbinned"
                host_status = "unresolved"
                _write_record(handle, contig_id, sequence)
            rows.append({"contig_id": contig_id, "genome_id": genome_id,
                         "binning_status": status, "host_assignment_status": host_status})
        if assignments:
            raise ValueError(f"Bin contigs absent from canonical assembly: {', '.join(sorted(assignments)[:5])}")
    return rows


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    """Write a table with a stable header, including for empty result sets."""
    with _atomic_text_output(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
