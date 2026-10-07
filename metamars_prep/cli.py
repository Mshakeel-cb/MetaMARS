"""The public metaMARS interface: sample-sheet inputs and resource controls."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

from .config import PreparationOptions
from .models import Sample, load_samples
from .taxonomy import RANKS


class HelpFormatter(argparse.HelpFormatter):
    """Keep flag names together; argument values appear in the usage line."""

    def __init__(self, prog: str):
        super().__init__(prog, max_help_position=26, width=110)

    def _format_action_invocation(self, action):
        return ", ".join(action.option_strings)


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        prog="metamars_prep", description="Prepare genomes and contigs for AMR evidence layers.",
        add_help=False, allow_abbrev=False, formatter_class=HelpFormatter,
    )
    command.add_argument("-h", "--help", action="help", help="help message")
    command.add_argument("-v", "--verbose", action="store_true", help="show commands and detailed tool output dflt=False")
    command.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    command.add_argument("-t", "--threads", type=positive_int, default=2, help="working threads dflt=2")
    command.add_argument("-j", "--jobs", type=positive_int, default=1, help="Maximum concurrent samples dflt=1")
    command.add_argument("-m", "--memory", type=positive_int, default=8, help="memory limit dflt=8")
    command.add_argument(
        "-tr", "--tax-rank", choices=RANKS, default="species", metavar="RANK",
        help="automatically group all classified genomes at specified rank "
             "{domain, phylum, class, order,family,genus,species} dflt=species",
    )
    command.add_argument(
        "-st", "--skip-taxonomy", action="store_true",
        help="Explicitly omit GTDB-Tk, groups are marked unresolved dflt=False",
    )
    command.add_argument("--db-dir", type=Path, default=Path.home() / "databases" / "metamars",
                         help="Directory for automatic database downloads dflt=~/databases/metamars")
    command.add_argument("--checkm2-db", help="Use an existing CheckM2 .dmnd file or directory")
    command.add_argument("--gtdbtk-db", help="Use an existing GTDB-Tk reference database directory")
    command.add_argument(
        "-it", "--input-type", metavar="TYPE",
        choices=("metagenome_reads", "isolate_reads", "isolate_genomes", "MAGs"),
        help="{metagenome_reads,isolate_reads,isolate_genomes,MAGs}",
    )
    command.add_argument(
        "-s", "--samples", type=Path, required=True,
        help="CSV/TSV sample sheet; paths are relative to this file",
    )
    return command


def input_samples(args) -> list[Sample]:
    # CheckM2's Prodigal wrapper and QUAST impose stricter path rules than Python.
    resolved_outdir = args.outdir.expanduser().resolve()
    if not re.fullmatch(r"[A-Za-z0-9_./-]+", str(resolved_outdir)):
        raise ValueError("--outdir must use letters, digits, '/', '.', '_' and '-' only; external tools do not safely support spaces or shell characters in output paths")
    return load_samples(args.samples.expanduser().resolve(), input_type=args.input_type)


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        samples = input_samples(args)
        from .pipeline import prepare
        prepare(samples, PreparationOptions(**vars(args)))
        return 0
    except KeyboardInterrupt:
        print("metaMARS interrupted. Logs and partial outputs are retained; use a new output directory to rerun.", file=sys.stderr)
        return 130
    except (ValueError, OSError, RuntimeError, EOFError) as exc:
        print(f"metaMARS error: {exc}", file=sys.stderr)
        return 1
