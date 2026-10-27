# Ground truth: the pure ODE (deng_cpg_ode, native simulate semantics)
# under the twin's exact perturbation (2 nA -> RG_HC_flx, 2.0-2.5 s) and
# the paper-exact 0.1 s window. Integer-index injection, 0.1 ms steps.
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
import deng_cpg_ode as D  # noqa: E402


def run_pert(pert_idx, t0, t1, dur=12.0, dt=0.1):
    n = int(dur * 1000 / dt)
    V = D.VREST.copy()
    h = np.ones(len(D.NAMES)) * (1 / (1 + D.KH * np.exp(D.SH * (D.EH - D.VREST_HC))))
    h[[12, 13]] = 0.0
    k0, k1 = int(100 / dt), int(102 / dt)      # 10 nA, 2 ms kickoff (twin's)
    p0, p1 = int(t0 * 1000 / dt), int(t1 * 1000 / dt)
    vE = np.zeros(n + 1); vF = np.zeros(n + 1)
    vE[0] = V[0]; vF[0] = V[1]
    for k in range(n):
        sat = np.clip((V[D.src] - D.elo) / (D.ehi - D.elo), 0, 1)
        isyn = np.zeros(len(D.NAMES))
        np.add.at(isyn, D.dst, D.g * sat * (D.esyn - V[D.dst]))
        m_inf = 1 / (1 + D.KM * np.exp(D.SM * (D.EM - V)))
        h_inf = 1 / (1 + D.KH * np.exp(D.SH * (D.EH - V)))
        i_na = np.zeros(len(D.NAMES))
        na = list(D.NA)
        i_na[na] = D.G_NA * m_inf[na] * h[na] ** D.KH * (D.E_NA - V[na])
        iapp = np.zeros(len(D.NAMES))
        if k0 <= k < k1:
            iapp[0] = 10.0
        if p0 <= k < p1:
            iapp[pert_idx] += 2.0
        V = V + (D.GM * (D.VREST - V) + iapp + isyn + i_na) * dt / D.CM
        h = h + dt * (h_inf - h) / D.TAU_H
        vE[k + 1] = V[0]; vF[k + 1] = V[1]
    return np.arange(n + 1) * dt / 1000.0, vE, vF


def onsets(t, v, thr=-45.0):
    on, last = [], -1e9
    above = v > thr
    for i in range(1, len(v)):
        if above[i] and not above[i - 1] and t[i] - last > 0.2 and t[i] > 1.0:
            on.append(round(float(t[i]), 3)); last = t[i]
    return on


FLX = D.NAMES.index("RG_HC_flx")
for tag, t0, t1 in (("0.1 s window", 2.0, 2.1), ("0.5 s window", 2.0, 2.5)):
    t, vE, vF = run_pert(FLX, t0, t1)
    oE = onsets(t, vE)
    inwin = [x for x in oE if t0 <= x <= t1 + 0.1]
    print(f"ODE 2 nA->RGflx {tag}: ext onsets {oE} | in-window {len(inwin)}"
          f" | flx max in-window {vF[(t>=t0)&(t<=t1)].max():.1f} mV")
