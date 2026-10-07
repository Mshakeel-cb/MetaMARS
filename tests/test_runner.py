"""Exercise subprocess execution with real tiny Python programs, without tools."""

import contextlib
import io
from pathlib import Path
import sys
import tempfile
import threading
import time
import unittest

from metamars_prep.runner import Runner, StageError


def python(source, *arguments):
    return [sys.executable, "-c", source, *arguments]


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.outdir = Path(self.temporary.name)
        self.cwd = self.outdir / "work" / "sample" / "stage"
        self.runner = Runner(self.outdir)
        self.addCleanup(self.runner.cancel)
        self.terminal = io.StringIO()
        redirect = contextlib.redirect_stderr(self.terminal)
        redirect.__enter__()
        self.addCleanup(redirect.__exit__, None, None, None)

    def run_stage(self, commands, outputs=()):
        return self.runner.run("sample.stage", commands, self.cwd, outputs=outputs)

    def test_stdout_stderr_streamed_and_saved(self):
        self.runner.verbose = True
        self.run_stage([python("import sys; from pathlib import Path; print('stdout message'); print('stderr message', file=sys.stderr); Path('result').write_text('data')")], outputs=["result"])
        self.assertEqual((self.cwd / "result").read_text(), "data")
        for message in ("stdout message", "stderr message"):
            self.assertIn(f"[sample/stage] {message}", self.terminal.getvalue())
            self.assertIn(message, (self.outdir / "logs/sample.stage.log").read_text())

    def test_pipeline_streams_data_and_logs_producer_stderr(self):
        self.run_stage([
            python("import sys; print('payload'); print('producer diagnostic', file=sys.stderr)"),
            python("import sys; from pathlib import Path; Path('result').write_text(sys.stdin.read())"),
        ], outputs=[self.cwd / "result"])
        self.assertEqual((self.cwd / "result").read_text(), "payload\n")
        self.assertIn("producer diagnostic", (self.outdir / "logs/sample.stage.log").read_text())

    def test_producer_failure_is_not_hidden_by_successful_consumer(self):
        with self.assertRaisesRegex(StageError, "command 1 exited 7"):
            self.run_stage([python("import sys; print('partial'); sys.exit(7)"), python("import sys; sys.stdin.read()")])

    def test_consumer_failure_stops_sleeping_producer(self):
        started = time.monotonic()
        with self.assertRaisesRegex(StageError, "command 2 exited 8"):
            self.run_stage([python("import time; time.sleep(60)"), python("import sys; sys.exit(8)")])
        self.assertLess(time.monotonic() - started, 10)

    def test_failed_consumer_startup_stops_producer_and_closes_pipes(self):
        started = time.monotonic()
        with self.assertRaises(FileNotFoundError):
            self.run_stage([python("import time; time.sleep(60)"), [str(self.outdir / "absent-executable")]])
        self.assertLess(time.monotonic() - started, 10)

    def test_missing_output_is_a_failure(self):
        with self.assertRaisesRegex(StageError, "did not produce"):
            self.run_stage([python("pass")], outputs=["missing"])

    def test_empty_bin_directory_is_valid(self):
        self.run_stage([python("from pathlib import Path; Path('bins').mkdir()")], outputs=["bins"])
        self.assertEqual(list((self.cwd / "bins").iterdir()), [])

    def test_arguments_are_not_shell_code_and_environment_is_passed(self):
        literal = "$(touch injected); `touch also_injected`"
        self.runner.run("sample.stage", [python("import os, sys; from pathlib import Path; Path('result').write_text(sys.argv[1] + os.environ['TEST_RUNNER_VALUE'])", literal)], self.cwd, outputs=["result"], env={"TEST_RUNNER_VALUE": "literal"})
        self.assertEqual((self.cwd / "result").read_text(), literal + "literal")
        self.assertFalse((self.cwd / "injected").exists())

    def test_commands_execute_every_time_without_checkpoints_or_cleanup(self):
        self.cwd.mkdir(parents=True)
        untouched = self.cwd / "unrelated"
        untouched.write_text("keep")
        command = python("from pathlib import Path; p = Path('count'); p.write_text(str(int(p.read_text()) + 1) if p.exists() else '1')")
        self.run_stage([command], outputs=["count"])
        self.run_stage([command], outputs=["count"])
        self.assertEqual((self.cwd / "count").read_text(), "2")
        self.assertEqual(untouched.read_text(), "keep")
        self.assertFalse((self.outdir / ".metamars").exists())

    def start_sleeper(self):
        errors = []
        marker = self.outdir / "sleeper-started"
        def execute():
            try:
                self.runner.run("other.stage", [python("from pathlib import Path; import sys,time; Path(sys.argv[1]).touch(); time.sleep(60)", str(marker))], self.outdir / "work" / "other")
            except StageError as error:
                errors.append(error)
        worker = threading.Thread(target=execute)
        worker.start()
        deadline = time.monotonic() + 5
        while not marker.exists() and time.monotonic() < deadline:
            time.sleep(0.02)
        self.assertTrue(marker.exists(), "Sleeper did not start")
        return worker, errors

    def test_cancel_stops_active_process_and_rejects_further_work(self):
        worker, errors = self.start_sleeper()
        self.runner.cancel()
        worker.join(timeout=10)
        self.assertFalse(worker.is_alive())
        self.assertTrue(errors)
        with self.assertRaisesRegex(StageError, "cancelled"):
            self.run_stage([python("pass")])

    def test_failure_cancels_another_concurrent_sample(self):
        worker, errors = self.start_sleeper()
        with self.assertRaisesRegex(StageError, "command 1 exited 9"):
            self.run_stage([python("import sys; sys.exit(9)")])
        worker.join(timeout=10)
        self.assertFalse(worker.is_alive())
        self.assertTrue(errors)


if __name__ == "__main__":
    unittest.main()
