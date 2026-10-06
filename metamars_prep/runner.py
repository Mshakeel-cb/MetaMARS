"""Run external commands with live output, stage logs, and cancellation."""

from __future__ import annotations

import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
import threading
import time


class StageError(RuntimeError):
    """A command failed or did not produce an expected output."""


class Runner:
    def __init__(self, outdir: Path) -> None:
        self.logs = Path(outdir).resolve() / "logs"
        self.logs.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        self._print_lock = threading.Lock()
        self._active: set[subprocess.Popen] = set()
        self._cancelled = threading.Event()

    def _show(self, stage: str, message: str) -> None:
        with self._print_lock:
            print(f"[{'/'.join(stage.rsplit('.', 1))}] {message}", file=sys.stderr, flush=True)

    def cancel(self) -> None:
        """Stop active commands and their process groups; reject further work."""
        self._cancelled.set()
        with self._lock:
            processes = list(self._active)
        for process in processes:
            try:
                if os.name == "posix":
                    os.killpg(process.pid, signal.SIGTERM)
                elif process.poll() is None:
                    process.terminate()
            except ProcessLookupError:
                pass
        deadline = time.monotonic() + 3
        for process in processes:
            try:
                process.wait(timeout=max(0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                pass
        for process in processes:
            try:
                if os.name == "posix":
                    # Descendants may survive after their group leader exits.
                    os.killpg(process.pid, signal.SIGKILL)
                elif process.poll() is None:
                    process.kill()
            except ProcessLookupError:
                pass
            process.wait()

    def run(self, stage: str, commands: list[list[str]], cwd: Path, *, outputs=(), env=None) -> None:
        """Execute argv lists as a pipeline without a shell.

        Stream every stderr and the last stdout; intermediate stdout is data
        for the next command. All exit codes and declared output paths are
        checked. A failure cancels every active stage. Nothing is cached.
        """
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,199}", stage):
            raise ValueError(f"Invalid stage name: {stage!r}")
        if not commands or any(not command for command in commands):
            raise ValueError("Provide at least one nonempty command")
        if self._cancelled.is_set():
            raise StageError("Run has been cancelled")
        cwd = Path(cwd).resolve()
        cwd.mkdir(parents=True, exist_ok=True)
        log_path = self.logs / f"{stage}.log"
        processes, readers, stream_errors = [], [], []
        log_lock = threading.Lock()
        started = time.monotonic()
        with log_path.open("w", encoding="utf-8") as log:
            def drain(stream, label):
                try:
                    for raw in iter(stream.readline, b""):
                        line = raw.decode("utf-8", errors="replace").rstrip("\n")
                        with log_lock:
                            log.write(f"[{label}] {line}\n")
                            log.flush()
                        self._show(stage, line)
                except Exception as error:
                    stream_errors.append(error)
                finally:
                    stream.close()

            try:
                rendered = " | ".join(shlex.join(command) for command in commands)
                log.write(f"Stage: {stage}\nWorking directory: {cwd}\nCommand: {rendered}\n")
                log.flush()
                self._show(stage, f"Running {rendered}")
                previous = None
                for index, command in enumerate(commands, 1):
                    with self._lock:
                        if self._cancelled.is_set():
                            raise StageError(f"{stage} was cancelled")
                        process = subprocess.Popen(
                            command, cwd=cwd, env={**os.environ, **(env or {})},
                            stdin=previous if previous is not None else subprocess.DEVNULL,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            start_new_session=os.name == "posix",
                        )
                        self._active.add(process)
                    processes.append(process)
                    if previous is not None:
                        previous.close()
                    previous = process.stdout
                    reader = threading.Thread(target=drain, args=(process.stderr, f"command {index} stderr"), daemon=True)
                    reader.start()
                    readers.append(reader)
                reader = threading.Thread(target=drain, args=(previous, "stdout"), daemon=True)
                reader.start()
                readers.append(reader)
                while any(process.poll() is None for process in processes):
                    if self._cancelled.is_set() or stream_errors or any(process.poll() not in (None, 0) for process in processes):
                        self.cancel()
                        break
                    time.sleep(0.05)
                for reader in readers:
                    reader.join()
                codes = [process.wait() for process in processes]
                failures = ", ".join(f"command {index} exited {code}" for index, code in enumerate(codes, 1) if code)
                if failures:
                    raise StageError(f"{stage} failed: {failures}. See {log_path}")
                if self._cancelled.is_set():
                    raise StageError(f"{stage} was cancelled. See {log_path}")
                if stream_errors:
                    raise StageError(f"Could not stream {stage} output: {stream_errors[0]}")
                for output in map(Path, outputs):
                    path = output if output.is_absolute() else cwd / output
                    if not path.exists():
                        raise StageError(f"{stage} did not produce {path}. See {log_path}")
            except BaseException:
                self.cancel()
                raise
            finally:
                for reader in readers:
                    reader.join()
                for process in processes:
                    for stream in (process.stdout, process.stderr):
                        if stream and not stream.closed:
                            stream.close()
                with self._lock:
                    self._active.difference_update(processes)
        self._show(stage, f"Finished in {time.monotonic() - started:.1f}s; log: {log_path}")
