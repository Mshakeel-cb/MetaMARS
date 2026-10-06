"""Public run settings and the preparation module's fixed scientific defaults."""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class PreparationOptions:
    samples: Path
    outdir: Path
    threads: int = 2
    jobs: int = 1
    memory: int = 8
    tax_rank: str = "species"
    skip_taxonomy: bool = False
    gtdbtk_db: str | None = None
    input_type: str | None = None

    min_completeness: float = field(init=False, default=50)
    max_contamination: float = field(init=False, default=10)
    min_bin_contig: int = field(init=False, default=2500)
    min_read_length: int = field(init=False, default=50)
    qualified_phred: int = field(init=False, default=20)
    min_read_pairs: int = field(init=False, default=1)
