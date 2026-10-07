"""Find tool commands on PATH or in the supplied Conda environments."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess

TOOL_NAMES = {
    "fastp": ("fastp",),
    "spades": ("spades.py", "spades"),
    "minimap2": ("minimap2",),
    "samtools": ("samtools",),
    "depth": ("jgi_summarize_bam_contig_depths",),
    "metabat2": ("metabat2",),
    "metaquast": ("metaquast.py", "metaquast"),
    "quast": ("quast.py", "quast"),
    "checkm2": ("checkm2",),
    "gtdbtk": ("gtdbtk",),
}

# These names match envs/preparation.yml, envs/checkm2.yml and envs/gtdbtk.yml.
TOOL_ENVIRONMENTS = {
    name: f"metamars-{name}" if name in {"checkm2", "gtdbtk"} else "metamars-preparation"
    for name in TOOL_NAMES
}


def required_tools(samples, skip_taxonomy: bool = False) -> set[str]:
    needed = {"checkm2"}
    if not skip_taxonomy:
        needed.add("gtdbtk")
    for sample in samples:
        if sample.input_type.endswith("_reads"):
            needed.update(("fastp", "spades", "minimap2", "samtools", "depth"))
        if sample.input_type == "metagenome_reads":
            needed.update(("metaquast", "metabat2"))
        else:
            needed.add("quast")
    return needed


def _find_executable(candidates: tuple[str, ...], path: str | None = None) -> str | None:
    for candidate in candidates:
        executable = shutil.which(candidate, path=path)
        if executable:
            return str(Path(executable).absolute())
    return None


def _environments(manager: str) -> list[Path]:
    """Ask an installed environment manager once; no activation or installation."""
    try:
        result = subprocess.run(
            [manager, "env", "list", "--json"], check=True,
            capture_output=True, text=True, timeout=20,
        )
        data = json.loads(result.stdout)
    except (OSError, subprocess.SubprocessError, ValueError):
        return []
    prefixes = data.get("envs", []) if isinstance(data, dict) else []
    if not isinstance(prefixes, list):
        return []
    return [Path(prefix).expanduser().absolute() for prefix in prefixes if isinstance(prefix, str) and prefix]


def discover_tools(names: set[str]) -> dict[str, list[str]]:
    """Return argv prefixes, preferring PATH over the named Conda environments.

    Tools in named environments run through the environment manager so their
    dependencies resolve correctly. Conda needs an explicit streaming flag;
    Mamba streams stdout/stderr by default.
    """
    unknown = names - TOOL_NAMES.keys()
    if unknown:
        raise ValueError("Unknown tool names: " + ", ".join(sorted(unknown)))
    commands = {}
    for name in sorted(names):
        executable = _find_executable(TOOL_NAMES[name])
        if executable:
            commands[name] = [executable]
    if len(commands) == len(names):
        return commands

    seen_managers = set()
    for candidate in ("conda", os.environ.get("CONDA_EXE"), "mamba", os.environ.get("MAMBA_EXE")):
        manager = _find_executable((candidate,)) if candidate else None
        if not manager or manager in seen_managers:
            continue
        seen_managers.add(manager)
        prefixes = _environments(manager)
        for name in sorted(names - commands.keys()):
            for prefix in prefixes:
                if prefix.name != TOOL_ENVIRONMENTS[name]:
                    continue
                executable = _find_executable(TOOL_NAMES[name], path=str(prefix / "bin"))
                if executable:
                    streaming = [] if "mamba" in Path(manager).name else ["--no-capture-output"]
                    commands[name] = [manager, "run", *streaming, "-p", str(prefix), executable]
                    break
        if len(commands) == len(names):
            return commands

    missing = sorted(names - commands.keys())
    environments = sorted({TOOL_ENVIRONMENTS[name] for name in missing})
    raise ValueError(
        "Missing tools: " + ", ".join(missing)
        + ". Add their executables to PATH or create the environments from envs/: "
        + ", ".join(environments) + ". Conda or Mamba must be available to run those environments."
    )
