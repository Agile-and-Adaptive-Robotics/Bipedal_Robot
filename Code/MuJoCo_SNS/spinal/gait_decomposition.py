"""Decomposed gait evaluation v3: uses the extracted winner jsons directly,
loads the curriculum stage params, and runs the runner for each model."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
MYO = r"C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
OUT_MD = SPINAL / "GAIT_DECOMPOSITION_20261006.md"
OUT_JSON = SPINAL / "GAIT_DECOMPOSITION_20261006.json"
JOINT_MAP = {"hip_flexion_r": "hip", "knee_angle_r": "knee", "ankle_angle_r": "ankle"}


def load_ref():
    sys.path.insert(0, str(SPINAL))
    from kine_ref import load_reference
    return load_reference()


def period_from_zc(x, dt):
    xc = np.where(np.diff(np.sign(x - np.mean(x))) != 0)[0]
    return float(np.mean(np.diff(xc)) * 2 * dt) if len(xc) >= 2 else None


def decompose(rj, rdt, sj, sdt):
    rj, sj = np.asarray(rj, float), np.asarray(sj, float)
    if len(rj) < 5 or len(sj) < 5:
        return {"error": "too short"}
    r_amp, s_amp = float(np.ptp(rj)), float(np.ptp(sj))
    r_per, s_per = period_from_zc(rj, rdt), period_from_zc(sj, sdt)
    n = max(len(rj), len(sj))
    r_n = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(rj)), rj)
    s_n = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(sj)), sj)
    r_n = (r_n - r_n.mean()) / (r_n.std() + 1e-12)
    s_n = (s_n - s_n.mean()) / (s_n.std() + 1e-12)
    xc = np.correlate(s_n, r_n, mode="full")
    peak = int(np.argmax(xc)) - (n - 1)
    corr = float(np.corrcoef(r_n, s_n)[0, 1])
    return {
        "amp_ratio": round(s_amp / r_amp, 3) if r_amp > 0.1 else None,
        "amp_robot": round(s_amp, 1),
        "amp_ref": round(r_amp, 1),
        "phase_lag_pct": round(peak / n * 100, 1),
        "corr_r": round(corr, 3),
        "freq_ratio": round(r_per / s_per, 3) if (r_per and s_per and s_per > 0) else None,
        "period_s": round(s_per, 3) if s_per else None,
    }


def main():
    ref = load_ref()
    ref_r = ref["r"]  # right-side reference: hip, knee, ankle arrays
    # Get ref dt from the kine_ref module
    # Reference is 100 points per cycle; T_r is the cycle period
    rdt = ref.get("T_r", 1.23) / 100.0

    MODELS = [
        ("s3k", "decomp_winner_s3k_mutfix.json"),
        ("syn6", "decomp_winner_syn6.json"),
        ("pruned", "decomp_winner_pruned.json"),
        ("w2lvar", "decomp_winner_w2lvar.json"),
    ]

    results = []
    for label, wj in MODELS:
        wpath = SPINAL / wj
        if not wpath.exists():
            results.append({"label": label, "error": f"{wj} not found"})
            continue
        w = json.loads(wpath.read_text())
        params = w.get("params", w)
        # Convert params to a flat dict the runner can consume
        flat = {}
        for k, v in params.items():
            if isinstance(v, dict):
                flat.update(v)
            else:
                flat[k] = v

        npz_out = SPINAL / f"_decomp_{label}.npz"
        # Build runner args: --fitted loads the fitted_walk_params base, then override
        # with the winner's stage params via --params-file
        args = [MYO, str(SPINAL / "runner.py"),
                "--fitted", "--eval", "--time", "22",
                "--params-file", str(wpath)]
        print(f"Running {label}...")
        try:
            import os
            env = {**os.environ, "RUNNER_DUMP": str(npz_out)}
            proc = subprocess.run(args, cwd=str(SPINAL), capture_output=True,
                                  text=True, timeout=600, encoding="utf8",
                                  errors="replace", env=env)
        except subprocess.TimeoutExpired:
            results.append({"label": label, "error": "runner timeout"})
            continue

        npz = npz_out if npz_out.exists() else SPINAL / "spinal_run.npz"
        if not npz.exists():
            results.append({"label": label, "error": f"no npz (rc={proc.returncode}): {proc.stdout[-200]}{proc.stderr[-200]}"})
            continue
        try:
            d = np.load(npz, allow_pickle=True)
        except Exception as e:
            results.append({"label": label, "error": f"npz load failed: {e}"})
            continue
        if "q" not in d or "key_joints" not in d:
            results.append({"label": label, "error": f"npz missing arrays; keys={list(d.keys())[:8]}"})
            continue

        q = d["q"]
        kj = [str(x) for x in d["key_joints"]]
        t = d["t"] if "t" in d else np.arange(len(q)) * 0.002
        sdt = float(np.median(np.diff(t)))

        entry = {"label": label, "score": w.get("score"), "joints": {}}
        for col, jname in JOINT_MAP.items():
            if col not in kj:
                left = col.replace("_r", "_l")
                if left in kj:
                    col = left
                else:
                    entry["joints"][jname] = {"error": "not in key_joints"}
                    continue
            sj = q[:, kj.index(col)]
            rj = ref_r.get(jname)
            if rj is None:
                entry["joints"][jname] = {"error": "no ref"}
                continue
            entry["joints"][jname] = decompose(rj, rdt, sj, sdt)
        results.append(entry)
        if npz != SPINAL / "spinal_run.npz":
            npz.unlink(missing_ok=True)

    OUT_JSON.write_text(json.dumps(results, indent=2, default=str), encoding="utf8")
    lines = ["# Gait decomposition — 2026-10-06", ""]
    lines.append("Score = kinematic error (less negative = better). "
                 "Amp ratio 1.0 = human-matched. Phase 0% = aligned. Freq ratio 1.0 = same cadence.")
    lines.append("")
    lines.append("| Model | Score | Joint | Amp ratio | Robot° | Ref° | Phase lag | Corr r | Freq ratio | Period |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in results:
        if "error" in r:
            lines.append(f"| {r['label']} | — | ALL | ERROR: {r['error'][:50]} | | | | | | |")
            continue
        score = f"{r['score']:.0f}" if isinstance(r.get("score"), (int, float)) else "—"
        for jname, jd in r.get("joints", {}).items():
            if "error" in jd:
                lines.append(f"| {r['label']} | {score} | {jname} | ERR | | | | | | |")
            else:
                lines.append(
                    f"| {r['label']} | {score} | {jname} "
                    f"| {jd.get('amp_ratio', '—')} "
                    f"| {jd.get('amp_robot', '—')} "
                    f"| {jd.get('amp_ref', '—')} "
                    f"| {jd.get('phase_lag_pct', '—')}% "
                    f"| {jd.get('corr_r', '—')} "
                    f"| {jd.get('freq_ratio', '—')} "
                    f"| {jd.get('period_s', '—')}s |"
                )
    OUT_MD.write_text("\n".join(lines), encoding="utf8")
    print(f"\nWrote {OUT_MD}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
