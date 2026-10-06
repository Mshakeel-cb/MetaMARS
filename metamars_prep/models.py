"""Validated input records shared by the CLI and preparation workflow."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field, replace
from pathlib import Path


INPUT_TYPES = frozenset({"metagenome_reads", "isolate_reads", "isolate_genomes", "MAGs"})
READ_TYPES = frozenset({"metagenome_reads", "isolate_reads"})
INPUT_FIELDS = frozenset({"sample_id", "input_type", "read1", "read2", "genome"})


@dataclass
class Sample:
    sample_id: str
    input_type: str
    read1: Path | None = None
    read2: Path | None = None
    genome: Path | None = None
    metadata: dict[str, str] = field(default_factory=dict)

    @property
    def files(self) -> list[Path]:
        return [path for path in (self.read1, self.read2, self.genome) if path is not None]


def validate_identifier(value: str, label: str = "sample_id") -> None:
    """Require IDs that can safely be used in FASTA headers and file names."""
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError(f"Invalid {label} {value!r}: use letters, digits, '.', '_' or '-', starting with a letter or digit")
    if len(value) > (128 if label == "sample_id" else 180):
        raise ValueError(f"{label} is too long for portable tool and output identifiers")


def _validate_file(value: Path | None, label: str) -> Path | None:
    if value is None:
        return None
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"{label} is not an existing regular file: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"{label} is empty: {path}")
    return path


def validate_sample(sample: Sample) -> Sample:
    """Return a validated copy with absolute paths; do not mutate the caller."""
    validate_identifier(sample.sample_id)
    if sample.input_type not in INPUT_TYPES:
        raise ValueError(f"Sample {sample.sample_id}: input_type must be one of {', '.join(sorted(INPUT_TYPES))}")
    if sample.input_type in READ_TYPES:
        if sample.read1 is None or sample.read2 is None:
            raise ValueError(f"Sample {sample.sample_id}: paired read1 and read2 are required for {sample.input_type}")
        if sample.genome is not None:
            raise ValueError(f"Sample {sample.sample_id}: genome cannot accompany read inputs")
    else:
        if sample.genome is None:
            raise ValueError(f"Sample {sample.sample_id}: genome is required for {sample.input_type}")
        if sample.read1 is not None or sample.read2 is not None:
            raise ValueError(f"Sample {sample.sample_id}: read inputs cannot accompany a supplied genome")
    result = replace(
        sample,
        read1=_validate_file(sample.read1, f"Sample {sample.sample_id} read1"),
        read2=_validate_file(sample.read2, f"Sample {sample.sample_id} read2"),
        genome=_validate_file(sample.genome, f"Sample {sample.sample_id} genome"),
        metadata=dict(sample.metadata),
    )
    if result.read1 is not None and result.read1.samefile(result.read2):
        raise ValueError(f"Sample {sample.sample_id}: read1 and read2 must be different files")
    return result


def load_samples(path: Path, input_type: str | None = None) -> list[Sample]:
    """Read a CSV/TSV manifest, optionally supplying a common input type.

    A supplied type fills missing or blank row types; explicit row types must
    agree with it. Without one, every row must declare its own type. Relative
    data paths belong to the manifest directory.
    """
    if input_type is not None and input_type not in INPUT_TYPES:
        raise ValueError(f"input_type must be one of {', '.join(sorted(INPUT_TYPES))}")
    manifest = Path(path).expanduser().resolve()
    if not manifest.is_file():
        raise ValueError(f"Input manifest is not an existing file: {manifest}")
    samples: list[Sample] = []
    seen: set[str] = set()
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        if manifest.suffix.lower() == ".csv":
            delimiter = ","
        elif manifest.suffix.lower() in {".tsv", ".tab"}:
            delimiter = "\t"
        else:
            preview = handle.read(8192)
            handle.seek(0)
            try:
                delimiter = csv.Sniffer().sniff(preview, delimiters=",\t").delimiter
            except csv.Error as error:
                raise ValueError("Input manifest must be a CSV or TSV table with a header") from error
        reader = csv.DictReader(handle, delimiter=delimiter, strict=True)
        try:
            fields = reader.fieldnames
        except csv.Error as error:
            raise ValueError(f"Malformed input manifest header: {error}") from error
        if not fields or "sample_id" not in fields:
            raise ValueError("Input manifest requires a sample_id column")
        if input_type is None and "input_type" not in fields:
            raise ValueError("Input manifest requires an input_type column when --input-type is not supplied")
        if any(not field.strip() for field in fields) or len(fields) != len(set(fields)):
            raise ValueError("Input manifest contains empty or duplicate column names")
        while True:
            try:
                row = next(reader)
            except StopIteration:
                break
            except csv.Error as error:
                raise ValueError(f"Manifest line {reader.line_num}: malformed CSV/TSV: {error}") from error
            if None in row:
                raise ValueError(f"Manifest line {reader.line_num}: more values than header columns")
            if not any(value and value.strip() for value in row.values()):
                continue
            values = {key: (value or "").strip() for key, value in row.items()}
            sample_id = values["sample_id"]
            if sample_id in seen:
                raise ValueError(f"Manifest line {reader.line_num}: duplicate sample_id {sample_id!r}")
            row_type = values.get("input_type", "")
            if input_type is not None and row_type and row_type != input_type:
                raise ValueError(
                    f"Manifest line {reader.line_num}: input_type {row_type!r} conflicts with --input-type {input_type!r}"
                )
            effective_type = row_type or input_type
            if effective_type is None:
                raise ValueError(f"Manifest line {reader.line_num}: input_type is required when --input-type is not supplied")

            def input_path(name: str) -> Path | None:
                value = values.get(name, "")
                if not value:
                    return None
                candidate = Path(value).expanduser()
                return candidate if candidate.is_absolute() else manifest.parent / candidate

            try:
                sample = validate_sample(Sample(
                    sample_id=sample_id,
                    input_type=effective_type,
                    read1=input_path("read1"),
                    read2=input_path("read2"),
                    genome=input_path("genome"),
                    metadata={key: value for key, value in values.items() if key not in INPUT_FIELDS},
                ))
            except ValueError as error:
                raise ValueError(f"Manifest line {reader.line_num}: {error}") from error
            samples.append(sample)
            seen.add(sample_id)
    if not samples:
        raise ValueError("Input manifest contains no samples")
    return samples
