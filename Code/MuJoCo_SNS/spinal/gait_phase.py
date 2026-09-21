"""Measured-GRF gait-event detection and stance/swing phase normalization.

The requested display convention maps each ipsilateral stance interval to
0--50% and the following swing interval to 50--100%, regardless of the
measured duty factor. This is intentionally a display normalization; the
original event times and duty factor remain available in the result.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[3]
OSIM_DIR = REPO / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
IK_MOT = OSIM_DIR / "subject01_walk1_ik.mot"
GRF_MOT = OSIM_DIR / "subject01_walk1_grf.mot"
PHASE = np.linspace(0.0, 100.0, 201)


def read_storage(path: Path):
    """Read an OpenSim .mot/.sto-like whitespace table."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    end = next(i for i, line in enumerate(lines)
               if line.strip().lower() == "endheader")
    names = lines[end + 1].split()
    rows = []
    for line in lines[end + 2:]:
        if not line.strip():
            continue
        try:
            rows.append([float(value) for value in line.split()])
        except ValueError:
            break
    values = np.asarray(rows, dtype=float)
    return values[:, 0], names[1:], values[:, 1:]


def contact_events(side: str, threshold_n: float = 50.0):
    """Return heel-strike and toe-off times from measured vertical GRF."""
    if side not in ("r", "l"):
        raise ValueError("side must be 'r' or 'l'")
    time, names, values = read_storage(GRF_MOT)
    prefix = "" if side == "r" else "1_"
    vertical = values[:, names.index(prefix + "ground_force_vy")]
    contact = vertical > threshold_n
    heel_strike = time[np.flatnonzero(np.diff(contact.astype(int)) == 1) + 1]
    toe_off = time[np.flatnonzero(np.diff(contact.astype(int)) == -1) + 1]
    return heel_strike, toe_off


def phase_normalize(time, values, side: str, grid=PHASE):
    """Normalize every complete measured stride in the supplied time range.

    Returns cycles [N, phase, channel], mean, standard deviation, measured
    duty factors, and the (heel strike, toe off, next heel strike) triples.
    Stance and swing each occupy half the output grid.
    """
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if len(time) != len(values):
        raise ValueError("time/value length mismatch")
    heel_strike, toe_off = contact_events(side)
    cycles, events, duties = [], [], []
    for start, stop in zip(heel_strike[:-1], heel_strike[1:]):
        if start < time[0] or stop > time[-1]:
            continue
        offs = toe_off[(toe_off > start) & (toe_off < stop)]
        if len(offs) != 1:
            continue
        off = float(offs[0])
        mask = (time >= start) & (time <= stop)
        if np.count_nonzero(mask) < 4:
            continue
        local_time = time[mask]
        phase = np.where(
            local_time <= off,
            50.0 * (local_time - start) / max(off - start, 1e-12),
            50.0 + 50.0 * (local_time - off) / max(stop - off, 1e-12),
        )
        cycles.append(np.column_stack([
            np.interp(grid, phase, values[mask, channel])
            for channel in range(values.shape[1])
        ]))
        events.append((float(start), off, float(stop)))
        duties.append((off - start) / (stop - start))
    if not cycles:
        raise RuntimeError(
            f"no complete {side}-side heel-strike/toe-off/heel-strike cycle "
            f"inside {time[0]:.3f}--{time[-1]:.3f} s")
    cycle_array = np.stack(cycles)
    return {
        "grid": np.asarray(grid),
        "cycles": cycle_array,
        "mean": cycle_array.mean(axis=0),
        "std": cycle_array.std(axis=0),
        "events": events,
        "duty": np.asarray(duties),
    }

