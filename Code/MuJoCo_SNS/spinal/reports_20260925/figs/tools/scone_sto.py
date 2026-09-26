"""Minimal SCONE/OpenSim .sto storage reader (header + whitespace table).

Shared parser for the 2026-09-25 walker-campaign figures. Read-only.
"""
import numpy as np


def read_sto(path):
    """Return (names, data) where data is an (nRows, nCols) float array.

    Header: lines up to 'endheader'; column-name line is the first line after
    it; data is whitespace/tab separated. Tolerates trailing blank lines.
    """
    with open(path, "r", errors="replace") as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip().lower() == "endheader":
            i += 1
            break
        i += 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    names = lines[i].split()
    i += 1
    rows = []
    for ln in lines[i:]:
        ln = ln.strip()
        if not ln:
            continue
        rows.append(ln.split())
    data = np.array(rows, dtype=float)
    if data.shape[1] != len(names):
        raise ValueError(
            "%s: %d names vs %d data columns" % (path, len(names), data.shape[1])
        )
    return names, data


def col(names, data, name):
    """Column by exact name; raises with the closest 5 candidates if missing."""
    if name in names:
        return data[:, names.index(name)]
    low = [n.lower() for n in names]
    if name.lower() in low:
        return data[:, low.index(name.lower())]
    import difflib

    cand = difflib.get_close_matches(name, names, n=5)
    raise KeyError("%s not found; closest: %s" % (name, cand))


if __name__ == "__main__":
    import sys

    for p in sys.argv[1:]:
        names, data = read_sto(p)
        t = col(names, data, "time")
        print("%s" % p)
        print("  rows=%d cols=%d t=[%.3f..%.3f]" % (data.shape[0], data.shape[1], t[0], t[-1]))
        for ax in ("pelvis.pos.x", "pelvis.pos.y", "pelvis.pos.z"):
            try:
                v = col(names, data, ax)
                print(
                    "  %-14s first=%8.3f last=%8.3f min=%8.3f max=%8.3f"
                    % (ax, v[0], v[-1], v.min(), v.max())
                )
            except KeyError as e:
                print("  %-14s MISSING (%s)" % (ax, e))
        print()
