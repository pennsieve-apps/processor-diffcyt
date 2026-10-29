"""The Lambda handler, without R: a stand-in Rscript on PATH records what it was given."""
import json
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from processor import handler  # noqa: E402

FAKE_RSCRIPT = """#!/bin/sh
python3 - "$@" <<'PY'
import json, os, sys
keep = {k: v for k, v in os.environ.items() if k in %s or k in ("INPUT_DIR", "OUTPUT_DIR", "TMPDIR")}
json.dump({"argv": sys.argv[1:], "env": keep}, open(os.environ["FAKE_RECORD"], "w"))
PY
exit ${FAKE_EXIT:-0}
"""


class HandlerTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        bin_dir = Path(self.tmp.name) / "bin"
        bin_dir.mkdir()
        script = bin_dir / "Rscript"
        script.write_text(FAKE_RSCRIPT % (list(handler.PARAMETERS),))
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
        self.record = Path(self.tmp.name) / "record.json"
        self.saved = dict(os.environ)
        os.environ["PATH"] = f"{bin_dir}{os.pathsep}{os.environ['PATH']}"
        os.environ["FAKE_RECORD"] = str(self.record)

    def tearDown(self):
        os.environ.clear()
        os.environ.update(self.saved)
        self.tmp.cleanup()

    def test_payload_becomes_environment(self):
        os.environ["DA_METHOD"] = "GLMM"  # left over from a previous warm invocation
        out = handler.handler({"inputDir": "/in", "outputDir": "/out", "FDR": 0.1, "MIN_CELLS": 5,
                               "BLOCK_COLUMN": "", "STATE_MARKERS": "CD25,HLADR"}, None)
        self.assertEqual(out, {"status": "ok", "outputDir": "/out"})
        rec = json.loads(self.record.read_text())
        self.assertTrue(rec["argv"][0].endswith("R/main.R"))
        env = rec["env"]
        self.assertEqual(env["INPUT_DIR"], "/in")
        self.assertEqual(env["OUTPUT_DIR"], "/out")
        self.assertEqual(env["FDR"], "0.1")
        self.assertEqual(env["MIN_CELLS"], "5")
        self.assertEqual(env["BLOCK_COLUMN"], "")  # empty = no block, kept as given
        self.assertEqual(env["STATE_MARKERS"], "CD25,HLADR")
        self.assertNotIn("DA_METHOD", env)  # not set this run -> R default
        self.assertEqual(env["TMPDIR"], "/tmp")

    def test_failure_raises(self):
        os.environ["FAKE_EXIT"] = "1"
        with self.assertRaises(RuntimeError):
            handler.handler({"inputDir": "/in", "outputDir": "/out"}, None)


if __name__ == "__main__":
    unittest.main()
