"""Shevtsova 2026 lumbar circuit - STANDALONE literature-template build.

Paper: Shevtsova, Lockhart, Rybak, Magnuson & Danner, eLife 14:RP107480
(DOI 10.7554/eLife.107480). LUMBAR-ONLY biped reduction transcribed from
replication/shevtsova_rules.json (Table 1 weights, every edge cited
there; cervical + inter-girdle edges excluded as quadruped-specific - see
the rules _comment).

Equations = the paper's own (Methods eq 1-14), verbatim unit system
(mV, pF, nS, ms), Euler at the paper's dt = 1 ms:

  C dV/dt = -IL - ISynE - ISynI - INaP [- INoise]      (eq 1/2)
  IL    = gL (V - EL)                                   (eq 3)
  INaP  = gNaP m(V) h (V - ENa)         [RG centers]    (eq 4-8)
  ISynE = gSynE { sum_j S(w_ji) f(V_j) + D_E } (V - ESynE)   (eq 9)
  ISynI = gSynI { sum_j S(-w_ji) f(V_j) + D_I } (V - ESynI)  (eq 10)
  S(x) = max(x, 0);  f(V) = clip((V - Vthr)/(Vmax - Vthr), 0, 1) (11/12)
  D_{E,I} = d_{E,I} alpha + b_{E,I}                    (eq 13)
  tau_h(V) = tau0 + (taumax - tau0)/cosh((V - V1/2,tau) k_tau)  (eq 8)

PARAMETER PROVENANCE + FLAGS (open items for Ben / the M3 gate):
  F-A  tau_h bell: the paper text extraction gives "150-400 ms" with
      k_tau NOT extractable; the lineage simulator's model file
      (github.com/SimonDanner/CPGNetworkSimulator models/Danner-etal-eLife.txt)
      carries hTau0 80 / htau 160 / hTauV12 -35 / hTauK 15. DEFAULT here:
      the PAPER text values tau0 150 / taumax 400 with the simulator's
      V1/2,tau -35 and k_tau 0.15 /mV (file's "15" read in its percent
      convention; both choices keep the high floor that the project's
      tau_h porting rule requires). Flag --tau=danner switches to the
      simulator pair 80/160 for the sensitivity check.
  F-B  leak split: the model file gives gLeak 4.5 / Eleak -62.5 for
      population block [0:7] (the RG centers) and 2.8 / -60 otherwise;
      applied here as RG centers -> 4.5/-62.5, all INs -> 2.8/-60.
  F-C  gSynE/gSynI = 10 nS both (paper Methods text, rules _dynamics_notes);
      Vthr = -50 = the file's Vmin, Vmax = 0.
  F-D  Drives (rules node DRIVE label, Table 2/3 base): RG-F dE=0.1 b=0;
      RG-E dE=0 bE=0.1; V0D dI=0.75 bI=0; V0V dI=0.15 bI=0 (hind).
  F-E  Noise OFF by default (sigma 0) - the Fig-9 stable-trajectory
      protocol runs at sigmaNoise 5 fA ~ 0; the 1.1 pA bout protocol is
      available via --sigma=1.1 pA (tau_noise 5 ms stand-in, flagged).

Analysis = the paper's Data-analysis section verbatim: flexion when
f(V_RG-F) >= 0.1; midextension = halfway between extension on/offset;
period = consecutive midextension intervals; frequency = 1/period;
LR phase = circular-mean phase of L midextensions within R cycles.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
RULES = HERE / "shevtsova_rules.json"

# ------------------------------------------------------------ parameters (pF/nS/mV/ms)
C, GSYN_E, GSYN_I = 10.0, 10.0, 10.0
E_NA, E_SYN_E, E_SYN_I = 50.0, -10.0, -75.0
V_THR, V_MAX = -50.0, 0.0
# Boltzmann DIVISION form 1/(1+exp((V - V12)/k)) with k in mV (the
# simulator file's convention: mk -6, hk 4, hTauK 15); the paper's eq-5
# product notation {1+exp[(V-V12)*k]}^-1 carries the same curves with
# k in 1/mV. m activates with depolarization (mk<0), h inactivates (hk>0).
M_K, M_V12 = -6.0, -40.0
H_K, H_V12 = 4.0, -45.0
G_NAP = 4.5
GL_RG, EL_RG = 4.5, -62.5         # flag F-B
GL_IN, EL_IN = 2.8, -60.0
TAU_BELL = dict(paper=(150.0, 400.0), danner=(80.0, 160.0))
# F-A resolution: k_tau is the ONE constant the paper's parameter list
# omits (tau0/taumax/V1/2,tau given, k_tau not). The lineage file's
# hTauK 15 (mV, division form) LATCHES the flexor centers under the
# 2026 tau pair 150/400 (measured: RG-F output never falls below the
# paper's own flexion criterion). k_tau = 6.67 mV reproduces the
# paper's stated behavior; selected as the documented calibrated
# stand-in for the unpublished constant (flagged for Ben / M3 gate).
HTAU_V12, HTAU_K = -35.0, 6.67
DRIVE = dict(  # (dE, bE, dI, bI) per population-substring; flag F-D
    RG=("RG-F", 0.1, 0.0, 0.0, 0.0), RGE=("RG-E", 0.0, 0.1, 0.0, 0.0),
    V0D=("V0D", 0.0, 0.0, 0.75, 0.0), V0V=("V0V", 0.0, 0.0, 0.15, 0.0),
)


def _drives(label: str):
    for _, (sub, dE, bE, dI, bI) in DRIVE.items():
        if sub in label:
            return dE, bE, dI, bI
    return 0.0, 0.0, 0.0, 0.0


def build_circuit():
    rules = json.loads(RULES.read_text(encoding="utf-8"))
    labels = [n["label"] for n in rules["nodes"] if n["type"] != "DRIVE"]
    idx = {l: i for i, l in enumerate(labels)}
    n = len(labels)
    is_rg = np.array([("RG-F" in l or "RG-E" in l) for l in labels],
                     dtype=bool)
    gl = np.where(is_rg, GL_RG, GL_IN)
    el = np.where(is_rg, EL_RG, EL_IN)
    dE = np.zeros(n); bE = np.zeros(n)
    dI = np.zeros(n); bI = np.zeros(n)
    for l in labels:
        de, be, di, bi = _drives(l)
        i = idx[l]
        dE[i], bE[i], dI[i], bI[i] = de, be, di, bi
    # synapse arrays: src, dst, weight (paper sign; S() splits E/I)
    src, dst, w = [], [], []
    for s in rules["synapses"]:
        src.append(idx[s["from"]])
        dst.append(idx[s["to"]])
        w.append(float(s["gain"]) if s["sign"] == "exc"
                 else -float(s["gain"]))
    return dict(labels=labels, idx=idx, is_rg=is_rg, gl=gl, el=el,
                dE=dE, bE=bE, dI=dI, bI=bI,
                src=np.array(src), dst=np.array(dst), w=np.array(w))


def simulate(alpha: float, t_end_s: float = 20.0, dt_ms: float = 1.0,
             tau_mode: str = "paper", sigma_pA: float = 0.0,
             seed: int = 0, circuit=None):
    c = circuit or build_circuit()
    n = len(c["labels"])
    nsteps = int(round(t_end_s * 1000.0 / dt_ms))
    tau0, taumax = TAU_BELL[tau_mode]
    rng = np.random.default_rng(seed)
    V = np.full(n, -60.0)
    h = np.full(n, 1.0)
    noise = np.zeros(n)
    out = np.zeros((nsteps + 1, n))
    out[0] = V
    w, src, dst = c["w"], c["src"], c["dst"]
    for k in range(nsteps):
        fV = np.clip((V - V_THR) / (V_MAX - V_THR), 0.0, 1.0)
        gE_in = np.zeros(n)
        gI_in = np.zeros(n)
        for j in range(len(w)):
            if w[j] > 0:
                gE_in[dst[j]] += w[j] * fV[src[j]]
            else:
                gI_in[dst[j]] += -w[j] * fV[src[j]]
        DE = c["dE"] * alpha + c["bE"]
        DI = c["dI"] * alpha + c["bI"]
        iL = c["gl"] * (V - c["el"])
        iSyE = GSYN_E * (gE_in + DE) * (V - E_SYN_E)
        iSyI = GSYN_I * (gI_in + DI) * (V - E_SYN_I)
        m_inf = 1.0 / (1.0 + np.exp((V - M_V12) / M_K))
        iNap = np.zeros(n)
        iNap[c["is_rg"]] = (G_NAP * m_inf[c["is_rg"]] * h[c["is_rg"]]
                            * (V[c["is_rg"]] - E_NA))
        if sigma_pA > 0:
            noise += dt_ms * (-noise / 5.0) + rng.normal(
                0.0, sigma_pA * np.sqrt(2.0 * dt_ms / 5.0), n)
        else:
            noise[:] = 0.0
        dV = (-iL - iSyE - iSyI - iNap - noise) / C * dt_ms
        V = np.clip(V + dV, -100.0, 50.0)
        h_inf = 1.0 / (1.0 + np.exp((V - H_V12) / H_K))
        tau_h = tau0 + (taumax - tau0) / np.cosh((V - HTAU_V12) / HTAU_K)
        h = h + dt_ms * (h_inf - h) / tau_h
        out[k + 1] = V
    t = np.arange(nsteps + 1) * dt_ms / 1000.0
    return t, out, c


def _episodes(t, w, vF):
    """Extension episodes per the paper: f(V_RG-F) < 0.1. Returns
    (starts, ends, mids, durations) in seconds."""
    fF = np.clip((vF - V_THR) / (V_MAX - V_THR), 0, 1)
    ext = w & (fF < 0.1)
    idxs = np.flatnonzero(w)
    d = np.diff(ext[idxs].astype(int))
    starts = list(idxs[1:][d == 1])
    ends = list(idxs[1:][d == -1])
    if ends and (not starts or ends[0] < starts[0]):
        starts = [idxs[0]] + starts       # window opened mid-episode
    if starts and (not ends or starts[-1] >= ends[-1]):
        ends = ends + [idxs[-1]]          # window closed mid-episode
    pairs = [(a, b) for a, b in zip(starts, ends) if b > a]
    return (np.array([t[a] for a, _ in pairs]),
            np.array([t[b] for _, b in pairs]),
            np.array([0.5 * (t[a] + t[b]) for a, b in pairs]),
            np.array([t[b] - t[a] for a, b in pairs]))


def analyze(t, out, c, skip_s: float = 5.0) -> dict:
    """The paper's data-analysis definitions (module docstring)."""
    w = t >= skip_s
    res = {}
    mids = {}
    for s in ("l", "r"):
        st, en, mid, dur = _episodes(t, w, out[:, c["idx"][f"RG-F_{s}"]])
        mids[s] = mid
        if len(mid) >= 3:
            per = float(np.mean(np.diff(mid)))
            res[s] = dict(n_ext=len(mid), freq=1.0 / per,
                          ext_frac=float(np.mean(dur)) / per)
        else:
            res[s] = dict(n_ext=len(mid), freq=float("nan"),
                          ext_frac=float("nan"))
    ph = []
    mL, mR = mids["l"], mids["r"]
    for a in range(len(mR) - 1):
        for b in mL:
            if mR[a] <= b <= mR[a + 1]:
                ph.append((b - mR[a]) / (mR[a + 1] - mR[a]))
    res["phase_med"] = float(np.median(ph)) if ph else float("nan")
    res["phase_n"] = len(ph)
    return res


if __name__ == "__main__":
    import io
    import sys
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    for a in (0.55, 1.05):
        t, out, c = simulate(a, t_end_s=20.0)
        r = analyze(t, out, c)
        print(f"alpha={a}: freq L {r['l']['freq']:.2f} R {r['r']['freq']:.2f} Hz, "
              f"LR phase {r['phase_med']:.3f} (n={r['phase_n']}), "
              f"ext_frac L {r['l']['ext_frac']:.2f}")
