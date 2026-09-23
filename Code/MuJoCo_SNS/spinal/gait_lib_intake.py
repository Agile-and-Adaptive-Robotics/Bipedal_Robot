"""Gait-library intake (2026-09-23): one command for Ben's SimTK downloads.

Drop the downloaded zips into  D:\\temp\\gait_lib_staging\\downloads\\
(default; pass another folder as argv[1]), then run this script (myo env):

    python gait_lib_intake.py [downloads_dir]

It (1) extracts every zip into a sibling folder, (2) walks the whole
staging tree for .mot files, classifying IK-style vs vertical-GRF files by
their column names, (3) pairs them per folder and builds a kine_ref-schema
reference via gait_lib_loader for every pair, saving each to
spinal/gait_refs/<name>.npz, and (4) writes a summary md listing what
integrated and what still needs manual pairing (e.g. muscfib
fiber-length files without GRF). Read-only w.r.t. core files.
"""
from __future__ import annotations

import io
import sys
import zipfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from gait_lib_loader import _read_mot, load_reference_general  # noqa: E402

STAGE = Path(r"D:\temp\gait_lib_staging")
REFS_OUT = HERE / "gait_refs"
IK_HINTS = ("hip_flexion", "knee_angle", "ankle_angle")
GRF_HINTS = ("ground_force_vy", "ground_force_vy")


def classify_mot(path: Path) -> str:
    try:
        _, names, _ = _read_mot(path)
    except Exception:
        return "unreadable"
    ns = " ".join(names)
    if any(h in ns for h in GRF_HINTS):
        return "grf"
    if any(h in ns for h in IK_HINTS):
        return "ik"
    return "other"


def main():
    dl = Path(sys.argv[1]) if len(sys.argv) > 1 else STAGE / "downloads"
    if dl.exists():
        zips = sorted(dl.glob("*.zip"))
        for z in zips:
            dst = dl / z.stem
            if not dst.exists():
                dst.mkdir(parents=True)
                with zipfile.ZipFile(z) as zf:
                    zf.extractall(dst)
                print(f"extracted {z.name} -> {dst}")
    mots = sorted(STAGE.rglob("*.mot"))
    print(f"{len(mots)} .mot files under {STAGE}")
    by_dir: dict[Path, dict[str, list[Path]]] = {}
    for m in mots:
        kind = classify_mot(m)
        by_dir.setdefault(m.parent, {"ik": [], "grf": [], "other": []})[
            kind].append(m)

    REFS_OUT.mkdir(exist_ok=True)
    lines = ["# Gait-library intake report (2026-09-23)", ""]
    n_ref = 0
    for d, kinds in sorted(by_dir.items()):
        pairs = []
        used = set()
        for ik in kinds["ik"]:
            for grf in kinds["grf"]:
                if grf in used:
                    continue
                pairs.append((ik, grf))
                used.add(grf)
                break
        for ik, grf in pairs:
            name = f"{d.name}_{ik.stem}"
            try:
                ref = load_reference_general(ik, grf)
            except Exception as e:
                lines.append(f"- FAILED {name}: {e}")
                continue
            out = {}
            for side in ("r", "l"):
                for j in ("hip", "knee", "ankle"):
                    out[f"{side}_{j}"] = ref[side][j]
            for k, v in ref.items():
                if k not in ("r", "l"):
                    out[k] = v
            np.savez(REFS_OUT / f"{name}.npz", **out)
            n_ref += 1
            lines.append(
                f"- OK {name}: T_r {ref['T_r']:.3f} s, duty "
                f"{ref['duty_r']:.2f}/{ref['duty_l']:.2f}, knee_min "
                f"{ref['knee_min_r']:.1f} deg, hip range "
                f"{ref['hip_range_r']:.1f} deg")
        other = kinds["other"] + [g for g in kinds["grf"]
                                  if g not in used]
        for o in other:
            lines.append(f"- UNPAIRED ({classify_mot(o)}): {o}")

    lines += ["", f"References saved: {n_ref} -> {REFS_OUT}",
              "Use: load a .npz and rebuild the dict (arrays r_/l_ + "
              "scalar keys) as kine_ref.REF_CACHE for scoring/training "
              "against it (see gait_lib_score.py / gait_lib_pilot2.py).",
              "", "muscfib fiber-length files (no GRF columns) will list "
              "as UNPAIRED - they integrate through the F-L-V validation "
              "route, not the reference route."]
    out_md = HERE / "reports_20260923" / "goal5_intake_report.md"
    out_md.parent.mkdir(exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"saved {out_md} ({n_ref} references)")


if __name__ == "__main__":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    main()
