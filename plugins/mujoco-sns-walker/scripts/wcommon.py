"""Shared helpers for the mujoco-sns-walker plugin scripts.

Machine-portable across the three AARL machines (see repo AGENTS.md machine
map). Everything here is stdlib-only so any python can import it; heavy
imports (optuna/numpy) stay in the scripts that need them, which are meant
to run under the SNS env (walker.cmd resolves it).

Resolution order everywhere: explicit env var override -> known per-machine
candidates -> sensible fallback. Never assume the plugin lives in the repo
(installed plugins are copied to a host cache).
"""
from __future__ import annotations

import io
import os
import socket
import subprocess
import sys
from pathlib import Path

REPO_CANDIDATES = [
    r"D:\Github\Bipedal_Robot",                       # EB475WS4
    r"D:\GitHub\Bipedal_Robot",                       # easteregg2 (case-insensitive FS makes both work)
    r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot",   # laptop DESKTOP-5Q16KE9
]

# The SNS/MuJoCo env on each machine (AGENTS.md machine map):
PY_CANDIDATES = [
    r"C:\Users\Ben Bolen\.conda\envs\myo\python.exe",    # EB475WS4 (env "myo")
    r"D:\Anaconda\envs\myo\python.exe",                  # easteregg2 (env "myo")
    r"C:\Users\Ben\.anaconda3\envs\myoconv\python.exe",  # laptop (env "myoconv")
]

DB_NAME = "optuna_walk.db"


def utf8_stdio() -> None:
    """Windows-console-safe stdout/stderr (copy-pasted ~95x in spinal/; here once)."""
    for stream in (sys.stdout, sys.stderr):
        try:
            if hasattr(stream, "buffer"):
                wrapped = io.TextIOWrapper(stream.buffer, encoding="utf-8", errors="replace")
                wrapped.reconfigure(line_buffering=True)  # type: ignore[attr-defined]
                if stream is sys.stdout:
                    sys.stdout = wrapped
                else:
                    sys.stderr = wrapped
        except Exception:
            pass


def resolve_spinal() -> Path:
    """Locate Code/MuJoCo_SNS/spinal: AARL_SPINAL env -> cwd -> repo candidates."""
    env = os.environ.get("AARL_SPINAL")
    if env:
        p = Path(env)
        if (p / "runner.py").exists():
            return p
    cwd = Path.cwd()
    if (cwd / "runner.py").exists():
        return cwd
    for repo in REPO_CANDIDATES:
        p = Path(repo) / "Code" / "MuJoCo_SNS" / "spinal"
        if (p / "runner.py").exists():
            return p
    raise SystemExit(
        "[walker] could not locate Code/MuJoCo_SNS/spinal - "
        "set AARL_SPINAL or run from the repo"
    )


def resolve_python() -> str:
    """Locate the SNS env python: AARL_PYTHON env -> per-machine candidates."""
    env = os.environ.get("AARL_PYTHON")
    if env and Path(env).exists():
        return env
    for cand in PY_CANDIDATES:
        if Path(cand).exists():
            return cand
    return sys.executable  # last resort: we may already BE the right env


def banner(msg: str = "") -> None:
    host = socket.gethostname()
    print(f"[walker @{host}] {msg}" if msg else f"[walker @{host}]")


def db_url(spinal: Path) -> str:
    """Optuna storage URL. Prefer running with cwd=spinal and the bare relative
    form (what the pipeline itself uses); absolute fallback for other cwds."""
    if Path.cwd().resolve() == spinal.resolve():
        return f"sqlite:///{DB_NAME}"
    return "sqlite:///" + (spinal / DB_NAME).as_posix()


MILESTONE_PATTERNS = (
    "== stage", "best ", "DONE", "FAILED", "seeded", "saved curriculum",
    "-> ", "applied",
)


def tail_milestones(log: Path, n: int = 3) -> list[str]:
    """Last n milestone-ish lines of a curriculum log (from the end, cheap)."""
    hits: list[str] = []
    try:
        with log.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                s = line.rstrip("\r\n")
                if any(p in s for p in MILESTONE_PATTERNS):
                    hits.append(s)
                    if len(hits) > 200:  # bound memory on big logs
                        hits.pop(0)
    except OSError:
        return []
    return hits[-n:]


def running_walker_processes() -> list[tuple[str, str, str]]:
    """(pid, name, commandline) of live processes that look like walker runs.
    Uses PowerShell CIM (tasklist has no command line)."""
    ps = (
        "$p = Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -match "
        "'runner\\.py|_curriculum\\.py|run_s3|run_curriculum|watch_walk' "
        "-and $_.CommandLine -notmatch 'Get-CimInstance' -and $_.ProcessId -ne $PID }; "
        "$p | ForEach-Object { \"$($_.ProcessId)`t$($_.Name)`t$($_.CommandLine)\" }"
    )
    try:
        out = subprocess.run(
            ["powershell", "-NoProfile", "-Command", ps],
            capture_output=True, text=True, timeout=30,
        ).stdout
    except Exception:
        return []
    procs = []
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3 and parts[0].strip().isdigit():
            procs.append((parts[0], parts[1], parts[2][:160]))
    return procs


def fmt_age(seconds: float) -> str:
    s = int(max(0, seconds))
    if s < 60:
        return f"{s}s ago"
    if s < 3600:
        return f"{s // 60}m ago"
    if s < 86400:
        return f"{s // 3600}h {(s % 3600) // 60}m ago"
    return f"{s // 86400}d ago"
