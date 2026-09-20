"""T1-T3: structured joint-layer PF basis fits (Ben, 2026-09-18).

T1 hypothesis: PF layers are JOINT/functional fields, each containing an
extensor and a flexor half-center:
    HIP   = {HIP-E, HIP-F}, KNEE = {KNEE-E, KNEE-F}, ANKLE = {ANKLE-E,
    ANKLE-F}   ->  6 HC waveforms, directly comparable to the
    unconstrained 6-synergy NMF.
  Group->HC mapping (anatomical): hip_ext/hip_abd -> HIP-E,
  hip_flex -> HIP-F, hip_add -> both (stance stabilizer vs pre-swing
  assist left to the fit), knee_ext -> KNEE-E, knee_flex -> KNEE-F,
  ankle_pf -> ANKLE-E, ankle_df -> ANKLE-F.
  Biarticular cross-terms: rect_fem +HIP-F, hamstrings +HIP-E,
  gastrocs +KNEE-F, gracilis +HIP-F, sartorius +KNEE-F.
T2: RF routing HIP-F vs KNEE-E vs both; sartorius likewise.
T3: merged layers: KNEE+ANKLE (4 HCs), HIP+KNEE (4 HCs), +TRUNK (8 HCs).

Protocol matches fsa_backsolve.py: 6 Hz lowpass + clip [0,1], muscle
filter (Fmax>5, max>0.05), sklearn NMF (nndsvda, 5000 it) as the
unconstrained reference, global-mean centered VAF, interleaved-frame
held-out (even=train, odd=test; test W solved by NNLS against the
FROZEN train H), 10 seeds.  Masked fits: masked multiplicative updates
(scale-stable) warm-started from the unconstrained solution.

Caveat (audit_ik_ground.py): activations back-solved on the
default-proportioned MJCF (~9 cm offset vs the scaled IK subject);
timings/co-activations usable, magnitudes carry the offset.

Outputs: stdout + fsa_results/joint_layers.json
"""
import io
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import nnls
from scipy.signal import butter, filtfilt
from sklearn.decomposition import NMF

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import gait_phase as gp
import muscle_map

OUT = Path("fsa_results/joint_layers.json")
SRC_NPZ = sys.argv[1] if len(sys.argv) > 1 else "bsolve_out.npz"
if SRC_NPZ != "bsolve_out.npz":
    _stem = Path(SRC_NPZ).stem.replace("bsolve_out", "").strip("_")
    OUT = Path(f"fsa_results/joint_layers_{_stem}.json")
LOWPASS_HZ = 6.0
N_SEEDS = 10
MU_ITERS = 4000

d = np.load(SRC_NPZ, allow_pickle=True)
acts = np.asarray(d["acts"], float)
names = [str(s) for s in d["act_names"]]
fmax = np.asarray(d["Fmax"], float)
t = np.asarray(d["t"], float)

dt = float(np.median(np.diff(t)))
bb, aa = butter(4, LOWPASS_HZ / (0.5 / dt))
SM = np.clip(filtfilt(bb, aa, acts, axis=0), 0.0, 1.0)

def centered_vaf(target, pred):
    sse = float(np.sum((target - pred) ** 2))
    sst = float(np.sum((target - target.mean()) ** 2))
    return 1.0 - sse / max(sst, 1e-12)

def phase_arr(side):
    hs_s, to_s = gp.contact_events(side)
    ph = np.full(len(t), np.nan)
    for h0, h1 in zip(hs_s[:-1], hs_s[1:]):
        m = (t >= h0) & (t <= h1)
        ph[m] = (t[m] - h0) / (h1 - h0) * 100.0
        for kk in to_s:
            if h0 < kk < h1:
                mm = (t >= h0) & (t <= kk)
                ph[mm] = (t[mm] - h0) / (kk - h0) * 50.0
                mm = (t > kk) & (t <= h1)
                ph[mm] = 50.0 + (t[mm] - kk) / (h1 - kk) * 50.0
    return ph

PH = {"r": phase_arr("r"), "l": phase_arr("l")}

HC6 = ("HIP-E", "HIP-F", "KNEE-E", "KNEE-F", "ANKLE-E", "ANKLE-F")
GROUP2HC = {"hip_ext": ("HIP-E",), "hip_abd": ("HIP-E",),
            "hip_add": ("HIP-E", "HIP-F"), "hip_flex": ("HIP-F",),
            "knee_ext": ("KNEE-E",), "knee_flex": ("KNEE-F",),
            "ankle_pf": ("ANKLE-E",), "ankle_df": ("ANKLE-F",)}
CROSS = {"rect_fem": ("HIP-F",), "semimem": ("HIP-E",),
         "semiten": ("HIP-E",), "bifemsh": ("HIP-E",),
         "gas_med": ("KNEE-F",), "gas_lat": ("KNEE-F",),
         "grac": ("HIP-F",), "sart": ("KNEE-F",)}

_cache = {}

def side_frames(side, with_trunk=False):
    key = (side, with_trunk)
    if key in _cache:
        return _cache[key]
    keep, lay = [], []
    for i, nm in enumerate(names):
        if not nm.endswith(f"_{side}") or fmax[i] <= 5.0 \
                or acts[:, i].max() <= 0.05:
            continue
        base = nm[:-2]
        g = muscle_map._GROUPS_BY_NAME.get(base)
        hcs = list(GROUP2HC.get(g[0], ())) if g else \
            (["TRUNK-E", "TRUNK-F"] if with_trunk else [])
        if not hcs:
            continue
        for sub, extra in CROSS.items():
            if base.startswith(sub):
                hcs += [x for x in extra if x not in hcs]
        keep.append(i)
        lay.append(hcs)
    r = (keep, lay, SM[:, keep])
    _cache[key] = r
    return r

def mask_of(keep, lay_allowed, hcs, overrides=None):
    M = np.zeros((len(hcs), len(lay_allowed)))
    for j, allowed in enumerate(lay_allowed):
        for L in allowed:
            if L in hcs:
                M[hcs.index(L), j] = 1.0
    for base, allowed in (overrides or {}).items():
        for j, gi in enumerate(keep):
            if names[gi][:-2].startswith(base):
                M[:, j] = 0.0
                for L in allowed:
                    if L in hcs:
                        M[hcs.index(L), j] = 1.0
    return M

def solve_w(A, H):
    W = np.zeros((A.shape[0], H.shape[0]))
    for i in range(A.shape[0]):
        W[i], _ = nnls(H.T, A[i])
    return W

def fit_masked_mu(A, M, H0, iters=MU_ITERS):
    """Masked multiplicative updates (scale-stable, monotone)."""
    H = np.maximum(H0, 1e-6) * M
    W = solve_w(A, H) + 1e-6
    eps = 1e-10
    for _ in range(iters):
        H *= (W.T @ A) / (W.T @ W @ H + eps)
        H *= M
        W *= (A @ H.T) / (W @ (H @ H.T) + eps)
    return H

def run_variant(side, hcs, overrides=None, label="", with_trunk=False,
                merge=None):
    keep, lay_allowed, A = side_frames(side, with_trunk=with_trunk)
    if merge:
        lay_allowed = [[merge.get(L, L) for L in al] for al in lay_allowed]
    M = mask_of(keep, lay_allowed, hcs, overrides)
    if (M.sum(axis=0) == 0).any():
        return {"error": "uncovered muscles"}
    k = len(hcs)
    tr = np.arange(0, A.shape[0], 2)
    te = np.arange(1, A.shape[0], 2)
    nmf = NMF(n_components=k, init="nndsvda", max_iter=5000, tol=1e-7,
              random_state=42)
    _Wfull = nmf.fit_transform(np.maximum(A, 0.0))
    Hfull = nmf.components_
    ves, best = [], None
    ph = PH[side]
    for s in range(N_SEEDS):
        rng = np.random.default_rng(500 + s)
        H0 = Hfull * (1.0 + 0.05 * rng.standard_normal(Hfull.shape))
        H = fit_masked_mu(A[tr], M, H0)
        Wte = solve_w(A[te], H)
        vaf = centered_vaf(A[te], Wte @ H)
        ves.append(float(vaf))
        if best is None or vaf > best[0]:
            Wf = solve_w(A, H)
            scale = np.maximum(Wf.max(axis=0), 1e-12)
            Wn = Wf / scale
            peaks = []
            for ci in range(k):
                cand = [j for j in range(A.shape[0]) if np.isfinite(ph[j])]
                jbest = cand[int(np.argmax(Wn[cand, ci]))]
                peaks.append(float(ph[jbest]))
            best = (vaf, Wn, H / scale.reshape(-1, 1), peaks)
    _, Wn, Hn, peaks = best
    return {"label": label, "side": side, "hcs": list(hcs),
            "muscles": [names[i] for i in keep],
            "vaf_held_mean": float(np.mean(ves)),
            "vaf_held_sd": float(np.std(ves)), "vaf_held_all": ves,
            "hc_peaks_pct": peaks, "H": Hn.tolist()}

results = {}
for side in ("r", "l"):
    keep, _, A = side_frames(side)
    tr, te = np.arange(0, A.shape[0], 2), np.arange(1, A.shape[0], 2)
    nmf = NMF(n_components=6, init="nndsvda", max_iter=5000, tol=1e-7,
              random_state=42)
    W6 = nmf.fit_transform(np.maximum(A, 0.0))
    H6 = nmf.components_
    Wte = solve_w(A[te], H6)
    results[f"FREE6_{side}"] = {
        "label": "unconstrained 6-comp NMF (reference)",
        "side": side, "vaf_in_sample": centered_vaf(A, W6 @ H6),
        "vaf_held_mean": centered_vaf(A[te], Wte @ H6), "vaf_held_sd": 0.0}

    results[f"T1_L3EF_{side}"] = run_variant(
        side, HC6, label="T1: 3 joint layers x E/F (6 HCs) + cross-terms")

    results[f"T3_KA_{side}"] = run_variant(
        side, ("HIP-E", "HIP-F", "KA-E", "KA-F"),
        merge={"KNEE-E": "KA-E", "KNEE-F": "KA-F",
               "ANKLE-E": "KA-E", "ANKLE-F": "KA-F"},
        overrides={"rect_fem": ("HIP-F", "KA-E"), "semimem": ("HIP-E",),
                   "semiten": ("HIP-E",), "bifemsh": ("HIP-E",),
                   "gas_med": ("KA-E", "KA-F"), "gas_lat": ("KA-E", "KA-F"),
                   "grac": ("HIP-F", "KA-F"), "sart": ("HIP-F", "KA-F")},
        label="T3: HIP + KNEE+ANKLE merged (4 HCs)")

    results[f"T3_HK_{side}"] = run_variant(
        side, ("HK-E", "HK-F", "ANKLE-E", "ANKLE-F"),
        merge={"HIP-E": "HK-E", "HIP-F": "HK-F",
               "KNEE-E": "HK-E", "KNEE-F": "HK-F"},
        overrides={"rect_fem": ("HK-E", "HK-F"),
                   "semimem": ("HK-E", "HK-F"),
                   "semiten": ("HK-E", "HK-F"),
                   "bifemsh": ("HK-E", "HK-F"),
                   "gas_med": ("HK-F", "ANKLE-E"),
                   "gas_lat": ("HK-F", "ANKLE-E"),
                   "grac": ("HK-F",), "sart": ("HK-F",)},
        label="T3: HIP+KNEE merged + ANKLE (4 HCs)")

    results[f"T3_TRUNK_{side}"] = run_variant(
        side, HC6 + ("TRUNK-E", "TRUNK-F"), with_trunk=True,
        label="T3: + TRUNK layer (8 HCs)")

    for tg, rf in (("rfHIPF", ("HIP-F",)), ("rfKNEEE", ("KNEE-E",)),
                   ("rfBOTH", ("KNEE-E", "HIP-F"))):
        results[f"{tg}_{side}"] = run_variant(
            side, HC6,
            overrides={"rect_fem": rf, "sart": ("KNEE-F", "HIP-F")},
            label=f"T2 routing {tag}" if False else f"T2 routing {tg}")

print("== Held-out centered VAF (interleaved frames, frozen-H, 10 seeds) ==")
for side in ("r", "l"):
    print(f"-- {side.upper()} side --")
    for tag in ("FREE6", "T1_L3EF", "T3_KA", "T3_HK", "T3_TRUNK",
                "rfHIPF", "rfKNEEE", "rfBOTH"):
        r = results.get(f"{tag}_{side}")
        if r and "vaf_held_mean" in r:
            extra = (f"  (in-sample {r['vaf_in_sample']:.3f})"
                     if "vaf_in_sample" in r else "")
            print(f"  {r['label']:<44s} {r['vaf_held_mean']:.3f} "
                  f"+- {r['vaf_held_sd']:.3f}{extra}")

for side in ("r", "l"):
    r = results[f"T1_L3EF_{side}"]
    print(f"\n== T1 {side.upper()}: HC peak phases (% cycle) ==")
    print("  " + ", ".join(f"{hc}={pk:.0f}%" for hc, pk in
                           zip(HC6, r["hc_peaks_pct"])))
    H = np.array(r["H"])
    print("  H column-norms: " + ", ".join(
        f"{hc}:{np.linalg.norm(H[:, j]):.2f}" for j, hc in enumerate(HC6)))

OUT.write_text(json.dumps(results, indent=1), encoding="utf-8")
print(f"\nwritten: {OUT}")
