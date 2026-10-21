"""Rybak et al. 2006a reduced two-level CPG - STANDALONE literature-template
build (ACTIVITY-BASED reinterpretation; a 1:1 Hodgkin-Huxley port is
explicitly NOT intended - replication/rybak_draft.md section 5).

Paper: Rybak, Shevtsova, Lafreniere-Roula & McCrea 2006a, "Modelling spinal
circuitry involved in locomotor pattern generation: insights from deletions
during fictive locomotion", J Physiol 577(Pt 2):617-639,
DOI 10.1113/jphysiol.2006.118703. Local full text:
spinal/lit_rybak2006_fulltext.txt (PMC1890439).

Circuit = the reduced unilateral hindlimb model of Fig 2A: 14 populations
(RG-E/F, Inrg-E/F, PF-E/F, Inpf-E/F, Ia-E/F, R-E/F, Mn-E/F) wired EXACTLY
per Table 2 through replication/rybak_rules.json (36 rows; every synapse
row in the json carries its cite - reproduced in the module CITES below).
NO invented weights: 0.0125, 0.45, 0.115, 0.0075, 0.05, 0.35, 0.2, 0.4,
0.1, 0.25, 0.5, 0.6, 0.2 are the paper's Table 2 values, verified verbatim
against the extracted table in the local full text.

Drives (d_mi; Table 2 weight "MLR (1)" x drive value, per figure legends):
  Fig 2B regime: d_rg-e = 0.5 > d_rg-f = 0.43  -> T_E > T_F
  Fig 2C regime: d_rg-e = 0.45 < d_rg-f = 0.51 -> T_F > T_E
  (quote, local text: "In B, the MLR drive to RG-E (d rg-e) population is
  larger than to RG-F population (d rg-f = 0.43; d rg-e = 0.5), and the
  model generates a rhythm with a longer duration extensor phase
  (T_E > T_F). In C, the RG-F population receives a larger drive
  (d rg-f = 0.51; d rg-e = 0.45), and the model generates a flexion-
  dominated rhythm (T_F > T_E).")
  d_pf = 0.5 to both PF populations ("These populations receive tonic
  excitatory drive from the MLR"; Table 2 row "PF-E: MLR (1)").
  Fig 6 resetting sims: d_rg-e = 0.5, d_rg-f = 0.46 (caption).
  Fig 10Ba non-resetting sims: d_rg-e = 0.55, d_rg-f = 0.44 (caption).

Equations - activity-based shell (the Rybak group's own later reduction of
this model family: Rybak 2015 eNeuro / Danner 2017 / Shevtsova 2026 - the
same shell as shevtsova_standalone.py in this folder), mV/pF/nS/ms:

  C dV/dt = -IL - ISynE - ISynI - INaP
  IL     = gL (V - EL)
  INaP   = gNaP m_inf(V) h (V - ENa)            [RG full; PF at paper ratio]
  ISynE  = (gSynE sum_j S(w_ji) f(V_j) + gDrive d_i + pulse) (V - ESynE)
  ISynI  =  gSynI sum_j S(-w_ji) f(V_j)         (V - ESynI)
  f(V)   = clip((V - Vthr)/(Vmax - Vthr), 0, 1);  S(x) = max(x, 0)

The paper's own eqn (10) separates the per-spike conductance increment
(gbar_E/gbar_I) from the DRIVE conductance (gbar_Ed/gbar_Id); the shell
keeps that separation (gSynE/gSynI vs gDrive).

INaP kinetics = paper Table 1 EXACT (extracted from the local full text
this session - supersedes rybak_draft.md's "tau_h,max NOT extracted"):
  m_inf = (1 + exp(-(V + 47.1)/3.1))^-1     instantaneous
  h_inf = (1 + exp( (V + 59)/8))^-1
  tau_h = 1200 / cosh((V + 59)/16)  ms      (tau_h,max = 1200 ms)
Division-form constants (Shevtsova-build lesson): km = -3.1 mV
(activation, NEGATIVE slope), kh = +8 mV (inactivation, POSITIVE slope).
tau_h floor in this build's operating range (plateaus approx -31..-15 mV
-> tau_h approx 154-405 ms) satisfies the project porting rule
(tau_h collapse quenches persistent-Na half-centers; DESIGN.md /
_tau_h_check.py).

FLAGS (documented stand-ins - the paper's appendix numeric table of
maximal conductances / leak / synaptic scalings did not survive the PMC
text conversion; "RG neurones: PF neurones: ..." fields extract empty):
  F-A gSynE = gSynI = GSYN = 90 nS - ONE scale multiplying ALL Table-2
      dimensionless weights; replaces the unextractable gbar_E/gbar_I
      x (20-neuron population rate) product of the spiking model.
      Calibrated to behavior per the rules json's own instruction
      ("Weights are ... NOT uS - scale by behavior when porting"):
      at 90 nS the paper's -0.115 Inrg inhibition holds a recovered RG
      center below escape (V_silent approx -65 mV), which is the
      paper's designed regime. NOT a weight change.
  F-B gDrive = 3.0 nS per unit d_mi - stand-in for the unextractable
      gbar_Ed; calibrated together with flag F-I so the paper's
      d = 0.43-0.55 keeps an uninhibited silent RG center just below
      the escape knee (approx -54 mV) and the inhibited center at
      approx -65 mV (m_inf < 0.01, full h recovery). Table 2 weights
      untouched.
  F-C Leak split: all non-RG populations 2.8 nS/-60 mV (lineage
      simulator convention, as in shevtsova_standalone.py flag F-B).
  F-I RG leak anchor EL_RG = -70 mV (GL_RG 4.5). Calibrated: the
      lineage file's -62.5 puts this paper's Table-1 INaP activation
      midpoint (-47.1 mV, k_m 3.1 mV) too close to the silent branch -
      the RG plateau then decays onto a stable middle fixed point
      (measured lock at V = -44.6, h = 0.14) instead of terminating.
      EL_RG = -70 (the paper's spiking-soma resting range for these
      exact Butera/Rybak-2003 INaP curves) places the inhibited silent
      branch at approx -65 mV where m_inf < 0.01 (true silence, h
      recovery to approx 0.68) and the burst terminates at the lower
      saddle-node knee, as in the paper's phasic regime.
  F-D gNaP: RG = 4.5 nS (calibrated burst/termination balance); PF =
      0.4 x RG - the paper's own ratio "0.1 mS cm-2 versus 0.25 mS
      cm-2" with PF "below the level allowing the expression of
      NaP-dependent rhythmogenic properties" (Methods) -> PF runs in
      forced mode and cannot oscillate alone, as in the paper.
  F-E Motoneurones: single-compartment activity-based populations; the
      Booth two-compartment structure and the MN-dendrite INaP are
      dropped per the draft's section-5 mapping (1:1 HH port not
      intended).
  F-G Table 2's Ia-F row prints "Ia-F (-0.1)" for the Ia-F target;
      wired as Ia-E <-> Ia-F MUTUAL inhibition per the rybak_rules.json
      cite (paper-typo flag, re-verified against the extracted table).
  F-H Deletion pulses are extra excitatory drive on the target
      population (paper's protocol: "increased excitatory drive"),
      shell-scaled in amplitude: PF-E pulse +1.0 d-unit (paper: "an
      88% increase in total drive", Fig 10Ba - their conductance
      regime differs), RG-E pulse +5.0 d-units (paper: "a temporary
      400% increase in excitatory drive to the RG-E population",
      Fig 6A). Behavioral criterion reproduced: PF pulse -> sustained
      PF-E + flexor MN deletion + rhythm phase maintained; RG pulse ->
      sustained RG-E ("changed RG-E population activity from phasic to
      sustained", Fig 6A) + phase shift. Pulse durations: 2.0 periods,
      inside the paper's 1-3-cycle deletion range (2 flexor bursts
      omitted, as in Fig 10Aa: "two bursts in the iTA trace were
      omitted").

Deletion-diagnostic quotes (local full text):
  Fig 6A: "the deletion of flexor activity was produced by a temporary
  400% increase in excitatory drive to the RG-E population ... The
  increased drive changed RG-E population activity from phasic to
  sustained ... when the applied 'perturbation' ends, the post-deletion
  rhythm re-appears with a phase shift in respect to the pre-deletion
  rhythm."
  Fig 8A/10Ba: "a deletion of flexor activity was produced by increased
  excitatory drive to the PF-E population ... the phase of the locomotor
  rhythm is maintained after the deletion."

Determinism: no noise, fixed ICs -> identical runs (the standing
"results repeat" directive); the verify script asserts trace identity.
Integration: explicit Euler at dt = 0.05 ms (< the fastest membrane
tau approx 0.083 ms at peak MN conductance).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
RULES = HERE / "rybak_rules.json"

# ------------------------------------------------ shell constants (flags F-A..F-E)
C = 10.0                       # pF
GSYN = 90.0                    # nS, single Table-2 weight scale (flag F-A)
GDRIVE = 3.0                   # nS per unit d_mi (flag F-B)
E_NA, E_SYN_E, E_SYN_I = 50.0, -10.0, -75.0
V_THR, V_MAX = -50.0, 0.0      # f(V) operating range (lineage shell)
G_NAP_RG = 4.5                 # nS (flag F-D)
G_NAP_PF = 0.4 * G_NAP_RG      # paper ratio 0.1/0.25 mS cm-2 (flag F-D)
GL_RG, EL_RG = 4.5, -70.0      # flags F-C/F-I (calibration note in docstring)
GL_IN, EL_IN = 2.8, -60.0

# --------------------------------- INaP kinetics: paper Table 1, division form
# m_inf = 1/(1+exp((V - M_V12)/M_K)), km NEGATIVE = activation (lesson 1)
M_K, M_V12 = -3.1, -47.1
# h_inf = 1/(1+exp((V - H_V12)/H_K)), kh POSITIVE = inactivation
H_K, H_V12 = 8.0, -59.0
# tau_h(V) = TAU_H_MAX / cosh((V - TAU_H_V12)/TAU_H_K), all paper values
TAU_H_MAX, TAU_H_V12, TAU_H_K = 1200.0, -59.0, 16.0

# Paper drive values (figure legends; see module docstring quotes)
DRIVE_REGIMES = {
    "fig2b_E_dominant": dict(d_rg_e=0.5, d_rg_f=0.43),
    "fig2c_F_dominant": dict(d_rg_e=0.45, d_rg_f=0.51),
    "fig6_resetting": dict(d_rg_e=0.5, d_rg_f=0.46),
    "fig10ba_nonresetting": dict(d_rg_e=0.55, d_rg_f=0.44),
}
D_PF = 0.5                     # MLR -> PF-E/PF-F tonic drive (Table 2 + text)


def build_circuit():
    """Wiring from rybak_rules.json (every row cited there).

    MLR rows are realized as the drive terms d_mi (Table 2 prints the
    weight "MLR (1)"; the drive VALUES come from the Fig 2/6/10 legends).
    """
    rules = json.loads(RULES.read_text(encoding="utf-8"))
    labels = [n["label"] for n in rules["nodes"] if n["type"] != "DRIVE"]
    idx = {l: i for i, l in enumerate(labels)}
    n = len(labels)
    is_rg = np.array(["RG-" in l for l in labels])
    is_pf = np.array(["PF-" in l for l in labels])
    gl = np.where(is_rg, GL_RG, GL_IN)
    el = np.where(is_rg, EL_RG, EL_IN)
    gnap = np.where(is_rg, G_NAP_RG, 0.0) + np.where(is_pf, G_NAP_PF, 0.0)
    WE = np.zeros((n, n))
    WI = np.zeros((n, n))
    cites = {}
    for s in rules["synapses"]:
        if s["from"] == "MLR":
            continue                    # drive rows enter as d_mi terms
        w = float(s["gain"]) * (1.0 if s["sign"] == "exc" else -1.0)
        if w > 0:
            WE[idx[s["to"]], idx[s["from"]]] += w
        else:
            WI[idx[s["to"]], idx[s["from"]]] += -w
        cites[(s["from"], s["to"])] = s["cite"]
    return dict(labels=labels, idx=idx, gl=gl, el=el, gnap=gnap,
                WE=WE, WI=WI, cites=cites, is_rg=is_rg)


def simulate(d_rg_e=0.5, d_rg_f=0.43, pulse=None, t_end_s=30.0,
             dt_ms=0.05, circuit=None):
    """Run the activity-based circuit. pulse = dict(target='PF-E'|'RG-E',
    t_start_s, dur_s, delta_d) adds extra excitatory drive (d-units x
    GDRIVE) to the target population during the window (flag F-H)."""
    c = circuit or build_circuit()
    n = len(c["labels"])
    nsteps = int(round(t_end_s * 1000.0 / dt_ms))
    d = np.zeros(n)
    d[c["idx"]["RG-E"]] = d_rg_e
    d[c["idx"]["RG-F"]] = d_rg_f
    d[c["idx"]["PF-E"]] = D_PF
    d[c["idx"]["PF-F"]] = D_PF
    gl, el, gnap = c["gl"], c["el"], c["gnap"]
    WE, WI = c["WE"], c["WI"]
    hasnap = gnap > 0
    pt = pulse or {}
    p_tgt = c["idx"][pt["target"]] if pt else -1
    p_t0 = pt.get("t_start_s", -1.0)
    p_t1 = p_t0 + pt.get("dur_s", 0.0)
    V = np.full(n, -60.0)
    h = np.full(n, 1.0)
    out = np.zeros((nsteps + 1, n))
    fout = np.zeros((nsteps + 1, n))
    out[0] = V
    for k in range(nsteps):
        t = k * dt_ms / 1000.0
        fV = np.clip((V - V_THR) / (V_MAX - V_THR), 0.0, 1.0)
        dd = d.copy()
        if pt and p_t0 <= t < p_t1:
            dd[p_tgt] += pt["delta_d"]
        gE_in = WE @ fV                # WE[to, from] -> row-vector times
        gI_in = WI @ fV
        iL = gl * (V - el)
        iSyE = (GSYN * gE_in + GDRIVE * dd) * (V - E_SYN_E)
        iSyI = GSYN * gI_in * (V - E_SYN_I)
        iNap = np.zeros(n)
        if hasnap.any():
            m_inf = 1.0 / (1.0 + np.exp((V - M_V12) / M_K))
            iNap[hasnap] = (gnap[hasnap] * m_inf[hasnap] * h[hasnap]
                            * (V[hasnap] - E_NA))
        V = np.clip(V + dt_ms * (-iL - iSyE - iSyI - iNap) / C, -100.0, 50.0)
        h_inf = 1.0 / (1.0 + np.exp((V - H_V12) / H_K))
        tau_h = TAU_H_MAX / np.cosh((V - TAU_H_V12) / TAU_H_K)
        h = h + dt_ms * (h_inf - h) / tau_h
        out[k + 1] = V
        fout[k + 1] = np.clip((V - V_THR) / (V_MAX - V_THR), 0.0, 1.0)
    t = np.arange(nsteps + 1) * dt_ms / 1000.0
    return t, out, fout, c


# ---------------------------------------------------- episode analysis (lesson 3)
def episodes(t, f, on=0.25, off=0.05):
    """Burst episodes of an f(V) trace with HYSTERETIC rising/falling edge
    pairing (a state machine pairs every rising edge with the NEXT falling
    edge; the last episode is marked partial if the window closes inside
    it). Returns list of dicts(t0, t1, partial)."""
    eps = []
    active = False
    t0 = None
    for i in range(1, len(f)):
        if not active and f[i] >= on and f[i - 1] < on:
            active, t0 = True, t[i]
        elif active and f[i] < off and f[i - 1] >= off:
            eps.append(dict(t0=t0, t1=t[i], partial=False))
            active = False
    if active:
        eps.append(dict(t0=t0, t1=t[-1], partial=True))
    return eps


def onsets(t, f, on=0.25, off=0.05):
    """Rising-edge onset times (s), same hysteresis as episodes()."""
    ts = []
    armed = True
    for i in range(1, len(f)):
        if armed and f[i] >= on and f[i - 1] < on:
            ts.append(t[i])
            armed = False
        elif not armed and f[i] < off:
            armed = True
    return np.array(ts)


def regime_stats(t, fout, c, skip_s=10.0):
    """Phase durations at the RG level (the rhythm's own phases) plus the
    MN-level pattern. T_E/T_F = mean RG-E/RG-F burst duration."""
    w = t >= skip_s
    res = {}
    for label in ("RG-E", "RG-F", "Mn-E", "Mn-F", "PF-E", "PF-F"):
        f = fout[:, c["idx"][label]]
        eps = [e for e in episodes(t, f) if e["t0"] >= skip_s
               and not e["partial"]]
        ons = [x for x in onsets(t, f) if x >= skip_s]
        durs = [e["t1"] - e["t0"] for e in eps]
        res[label] = dict(n=len(eps), dur_mean=float(np.mean(durs)) if durs
                          else float("nan"),
                          dur_all=durs, onsets=np.array(ons))
    rge = res["RG-E"]["onsets"]
    res["period"] = (float(np.mean(np.diff(rge))) if len(rge) >= 3
                     else float("nan"))
    res["period_cv"] = (float(np.std(np.diff(rge)) / np.mean(np.diff(rge)))
                        if len(rge) >= 3 else float("nan"))
    return res


if __name__ == "__main__":
    import io
    import sys
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    for name in ("fig2b_E_dominant", "fig2c_F_dominant"):
        dr = DRIVE_REGIMES[name]
        t, out, fout, c = simulate(**dr, t_end_s=30.0)
        st = regime_stats(t, fout, c)
        print(f"{name} (d_rg-e={dr['d_rg_e']}, d_rg-f={dr['d_rg_f']}): "
              f"period {st['period']:.3f} s (cv {st['period_cv']:.3f}), "
              f"T_E {st['RG-E']['dur_mean']:.3f} s, "
              f"T_F {st['RG-F']['dur_mean']:.3f} s, "
              f"MN-E {st['Mn-E']['dur_mean']:.3f} s, "
              f"MN-F {st['Mn-F']['dur_mean']:.3f} s "
              f"[{st['RG-E']['n']}/{st['RG-F']['n']} RG eps]")
