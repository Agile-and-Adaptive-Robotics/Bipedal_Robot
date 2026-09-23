"""Launch gate for walker runs: compile check + env-resolved launcher.

The session pattern (cd + set CONDA_PREFIX + full-path python + >> log) is
fragile through agent shells - compound cmd chains silently fail. This script
does the resolution in Python and either PRINTS the exact launch info or
launches detached itself with env vars set in-process (no shell chaining).

Usage:
  run_gate.py compile                          # py_compile the 4 core files
  run_gate.py print runner --eval --drive 2.5  # show resolved launch info
  run_gate.py print curriculum 3 40
  run_gate.py go <tag> runner --eval --drive 2.5     # detached, log <tag>_YYYYMMDD.log
  run_gate.py go <tag> curriculum 3 40
Notes:
  - `go` writes to spinal/<tag>_<YYYYMMDD>.log and prints PID + log path.
  - Unless the runner args contain an AARL_NPZ-worthy name already, `go`
    sets AARL_NPZ=spinal_run_<tag>.npz (unique name -> no Spyder/AV file locks
    on spinal_run.npz).
"""
from __future__ import annotations

import os
import py_compile
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import wcommon  # noqa: E402

CORE_FILES = ("runner.py", "_curriculum.py", "params.py", "build_network.py")


def compile_gate(spinal: Path) -> int:
    wcommon.banner(f"compile gate in {spinal}")
    bad = 0
    for f in CORE_FILES:
        try:
            py_compile.compile(str(spinal / f), doraise=True)
            print(f"  ok    {f}")
        except py_compile.PyCompileError as e:
            bad += 1
            print(f"  FAIL  {f}: {e}")
    print("compile gate: PASS" if not bad else "compile gate: FAIL - do not launch")
    return 1 if bad else 0


def build_launch(spinal: Path, mode: str, rest: list[str], tag: str | None):
    py = wcommon.resolve_python()
    env = {**os.environ, "CONDA_PREFIX": str(Path(py).parent)}
    npz = None
    if mode == "runner":
        argv = [py, "runner.py"] + rest
        if tag:
            npz = f"spinal_run_{tag}.npz"
    elif mode == "curriculum":
        stage = rest[0] if rest else "3"
        n = rest[1] if len(rest) > 1 else "40"
        argv = [py, "_curriculum.py", stage, n]
    else:
        raise SystemExit(f"unknown mode '{mode}' (runner|curriculum)")
    if npz:
        env["AARL_NPZ"] = npz
    return py, argv, env, npz


def main() -> int:
    wcommon.utf8_stdio()
    # Manual prefix parse: run_gate.py <action> [tag] <mode> [passthrough args...]
    # (argparse would eat runner flags like --eval / --drive).
    argv = sys.argv[1:]
    if not argv or argv[0] not in ("compile", "print", "go"):
        print(__doc__)
        return 0 if not argv else 1
    action = argv[0]
    rest = argv[1:]
    tag = None
    if action in ("print", "go") and rest and rest[0] not in ("runner", "curriculum"):
        tag = rest.pop(0)
    if not rest or rest[0] not in ("runner", "curriculum"):
        print(f"{action} needs a mode: runner|curriculum", file=sys.stderr)
        return 1
    mode = rest.pop(0)

    spinal = wcommon.resolve_spinal()

    if action == "compile":
        return compile_gate(spinal)

    py, argv_run, env, npz = build_launch(spinal, mode, rest, tag)
    wcommon.banner(f"spinal={spinal}")
    print(f"python : {py}")
    print(f"cwd    : {spinal}")
    print(f"argv   : {' '.join(argv_run[1:])}")
    print(f"env    : CONDA_PREFIX={env['CONDA_PREFIX']}"
          + (f" AARL_NPZ={npz}" if npz else ""))

    if action == "print":
        print("(dry run - nothing launched)")
        return 0

    if not tag:
        print("go needs a tag (used for log + npz names)", file=sys.stderr)
        return 1
    log = spinal / f"{tag}_{time.strftime('%Y%m%d')}.log"
    print(f"log    : {log}")
    creationflags = 0
    if os.name == "nt":
        creationflags = subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP
    with log.open("ab") as fh:
        fh.write(f"\n=== {tag} start {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n".encode())
        fh.write(f"=== argv: {argv_run[1:]} ===\n".encode())
        fh.flush()
        proc = subprocess.Popen(
            argv_run, cwd=str(spinal), env=env,
            stdout=fh, stderr=subprocess.STDOUT,
            creationflags=creationflags,
        )
    print(f"launched detached: PID {proc.pid}")
    print("poll with: walker.cmd status.py   (or tail the log)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
