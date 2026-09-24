"""Gait-library intake v2 (2026-09-24): thumb-drive edition.

Handles the F:\\Biomechanics data haul staged into
D:\\temp\\gait_lib_staging\\downloads\\:
  - nested zips (SubjectNN-latest.zip -> subjectNN.zip -> tree; the
    inner zip is extracted recursively, __MACOSX skipped)
  - Arnold-style layouts where IK lives in ik/Results_191/Run_XXXXX.mot
    and GRF in ExportedData/Run_X XX_newCOP3.mot (cross-FOLDER pairing
    by trial number)
  - plain flat layouts (Falisse Case_40: per-folder pairing, as v1)

Run (myo env): python gait_lib_intake.py [downloads_dir]
Writes refs to spinal/gait_refs/<name>.npz (kine_ref schema) and the
summary md to reports_20260923/goal5_intake_report.md.
Read-only w.r.t. core files; staging stays outside the repo.
"""
from __future__ import annotations

import io
import re
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
GRF_HINTS = ("ground_force",)          # widened: plate exports vary


def extract_recursive(root: Path, depth: int = 0) -> None:
    if depth > 4:
        return
    for z in sorted(root.rglob("*.zip")):
        if "__MACOSX" in str(z):
            continue
        dst = z.with_suffix("")          # subjectNN.zip -> subjectNN/
        if dst.exists():
            continue
        try:
            with zipfile.ZipFile(z) as zf:
                members = [m for m in zf.namelist()
                           if "__MACOSX" not in m]
                zf.extractall(dst, members=members)
            print(f"extracted {z.name} -> {dst}")
            extract_recursive(dst, depth + 1)
        except zipfile.BadZipFile:
            print(f"BAD ZIP skipped: {z}")


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


def trial_token(name: str) -> str | None:
    """Run_20002 / Run_200 02_newCOP3 / Run_30002_v24 -> '20002'."""
    m = re.search(r"[Rr]un[_ ]?(\d{2,5})\s?(\d{2})?", name)
    if not m:
        return None
    return (m.group(1) + (m.group(2) or "")).lstrip("0") or "0"


def main():
    dl = Path(sys.argv[1]) if len(sys.argv) > 1 else STAGE / "downloads"
    if dl.exists():
        extract_recursive(dl)
    mots = [m for m in sorted(STAGE.rglob("*.mot"))
            if "__MACOSX" not in str(m)]
    print(f"{len(mots)} .mot files under {STAGE}")
    by_dir: dict[Path, dict[str, list[Path]]] = {}
    for m in mots:
        kind = classify_mot(m)
        if kind == "unreadable":
            kind = "other"
        by_dir.setdefault(m.parent, {"ik": [], "grf": [], "other": []})[
            kind].append(m)

    REFS_OUT.mkdir(exist_ok=True)
    lines = ["# Gait-library intake report (2026-09-24, thumb drive)", ""]
    n_ref = 0

    def build_ref(name: str, ik: Path, grf: Path) -> None:
        nonlocal n_ref
        try:
            ref = load_reference_general(ik, grf)
        except Exception as e:
            lines.append(f"- FAILED {name}: {e}")
            return
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

    # ---- pass 1: per-folder pairing (Falisse-style flat layouts) ----
    # Skip RRA/CMC derivative trees: their states files carry angles AND
    # GRF columns, so they self-pair into duplicates of the clean IK.
    def is_derivative(p: Path) -> bool:
        s = str(p).lower()
        return ("rra" in s.split("\\") or "cmc" in s.split("/")
                or "/rra" in s or "/cmc" in s or "\\rra" in s
                or "\\cmc" in s)

    handled: set[Path] = set()
    for d, kinds in sorted(by_dir.items()):
        if is_derivative(d):
            continue
        pairs, used = [], set()
        for ik in kinds["ik"]:
            for grf in kinds["grf"]:
                if grf in used:
                    continue
                pairs.append((ik, grf))
                used.add(grf)
                break
        for ik, grf in pairs:
            build_ref(f"{d.name}_{ik.stem}", ik, grf)
            handled.update({ik, grf})

    # ---- pass 2: cross-folder trial pairing inside one subject tree ----
    # group folders by their nearest "subjectNN" ancestor (or zip stem)
    def subject_root(p: Path) -> Path | None:
        for part in p.parents:
            if re.match(r"subject\d+", part.name, re.I):
                return part
            if part == STAGE:
                return None
        return None

    groups: dict[Path, dict[str, list[Path]]] = {}
    for d, kinds in by_dir.items():
        root = subject_root(d)
        if root is None:
            continue
        g = groups.setdefault(root, {"ik": [], "grf": []})
        for kind in ("ik", "grf"):
            for f in kinds[kind]:
                if f not in handled:
                    g[kind].append(f)

    for root, g in sorted(groups.items()):
        # keep clean IK sources (ik/Results_191/*.mot); drop RRA/CMC
        # per-cycle derivatives (Kinematics_q / states_degrees)
        clean_ik = [f for f in g["ik"]
                    if "kinematics_q" not in f.name.lower()
                    and "states_degrees" not in f.name.lower()
                    and not is_derivative(f)]
        pool_ik = clean_ik or [f for f in g["ik"]
                               if "kinematics_q" not in f.name.lower()
                               and "states_degrees" not in f.name.lower()]
        grf_by_tok = {}
        for grf in g["grf"]:
            tok = trial_token(grf.name)
            if tok is None:
                continue
            # prefer non-_v24 exports when both exist for a token
            if tok not in grf_by_tok or (
                    "_v24" in grf_by_tok[tok].name and "_v24" not in grf.name):
                grf_by_tok[tok] = grf
        for ik in pool_ik:
            tok = trial_token(ik.name)
            grf = grf_by_tok.get(tok) if tok else None
            if grf is None and len(pool_ik) == 1 and g["grf"]:
                # single IK in the group (e.g. running): pick the GRF
                # with the most name-word overlap
                words = lambda p: set(re.findall(r"[a-z]+", p.name.lower()))
                best = max(g["grf"], key=lambda gp: len(words(gp) & words(ik)))
                if len(words(best) & words(ik)) >= 1:
                    grf = best
            if grf is None:
                continue
            build_ref(f"{root.name}_{ik.stem}", ik, grf)
            handled.update({ik, grf})

    # ---- leftovers report ----
    n_unpaired = 0
    for d, kinds in sorted(by_dir.items()):
        other = [f for f in kinds["other"] + kinds["ik"] + kinds["grf"]
                 if f not in handled]
        for o in other[:6]:                    # cap the spam per folder
            lines.append(f"- UNPAIRED ({classify_mot(o)}): {o}")
            n_unpaired += 1
        if len(other) > 6:
            lines.append(f"  ... +{len(other) - 6} more in {d}")

    lines += ["", f"References saved: {n_ref} -> {REFS_OUT}",
              f"Unpaired/leftover .mot listed: {n_unpaired}",
              "Use: load a .npz and rebuild the dict (arrays r_/l_ + "
              "scalar keys) as kine_ref.REF_CACHE for scoring/training "
              "against it (see gait_lib_score.py / gait_lib_pilot2.py).",
              "", "Predictive .sto libraries (results-speeds / deficits / "
              "assistloadwalk states) integrate through a separate "
              "states-based route, not this GRF-phased one."]
    out_md = HERE / "reports_20260923" / "goal5_intake_report.md"
    out_md.parent.mkdir(exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"saved {out_md} ({n_ref} references, {n_unpaired} unpaired)")


if __name__ == "__main__":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    main()
