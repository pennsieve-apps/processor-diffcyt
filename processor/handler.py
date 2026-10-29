"""
The Lambda runtime. A lambda step's payload carries inputDir, outputDir and
the step's parameters as fields (where a container gets them as environment
variables); the handler sets those variables and runs the same R script as a
subprocess, so both runtimes write the same outputs.
"""
import json
import os
import subprocess
from pathlib import Path

PARAMETERS = (
    "CONDITION_COLUMN",
    "REFERENCE_LEVEL",
    "BLOCK_COLUMN",
    "DA_METHOD",
    "DS_METHOD",
    "MIN_CELLS",
    "MIN_SAMPLES",
    "FDR",
    "STATE_MARKERS",
)


def main_script() -> str:
    root = os.environ.get("APP_ROOT") or str(Path(__file__).resolve().parent.parent)
    return str(Path(root) / "R" / "main.R")


def _env_value(value) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (dict, list)):
        return json.dumps(value)
    return str(value)


def handler(event, context):
    env = dict(os.environ)
    env["INPUT_DIR"] = event.get("inputDir") or ""
    env["OUTPUT_DIR"] = event.get("outputDir") or ""
    # A warm function keeps its environment between invocations; building the
    # child's environment per call means a parameter this run doesn't set falls
    # back to its default instead of the last run's value.
    for name in PARAMETERS:
        env.pop(name, None)
        value = event.get(name)
        if value is not None:
            env[name] = _env_value(value)
    # R and its packages write nothing outside /tmp, the only writable path on Lambda.
    env.setdefault("HOME", "/tmp")
    env["TMPDIR"] = "/tmp"
    # Output streams straight to the function's log (CloudWatch).
    proc = subprocess.run(["Rscript", main_script()], env=env, check=False)
    if proc.returncode != 0:
        # An exception, so the step fails; a returned value would succeed.
        raise RuntimeError(f"processor-diffcyt exited with {proc.returncode}")
    return {"status": "ok", "outputDir": env["OUTPUT_DIR"]}
