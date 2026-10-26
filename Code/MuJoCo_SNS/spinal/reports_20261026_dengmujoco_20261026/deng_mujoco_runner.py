"""Deng 2023 perturbation protocol, MuJoCo replication (EB475WS4, 2026-10-26).

Fourth leg of the comparison: paper / AnimatLab original (easteregg2,
2026-10-26) / Simulink twin (this machine) / THIS runner.

Faithfulness (per prompt ruling): the Deng organization has NO contact
sensors and the 2023 paper suspends the pelvis in air - so this is the
Nourse-template circuit driving the W2L body with the root WELDED
(w2l_air.xml: freejoint removed = constant-height hold), LEFT leg only
(the rat model is one hindlimb; right leg hangs inert, ctrl 0), NO
contact feedback, and NO afferent glue - deng_cpg_ode carries no Ia/II/Ib
populations, so the coupling is feedforward, which is exactly the
template's content. Documented substitution: the template has MNs only
for hip; knee/ankle are driven from the PF_ka half-centers through the
same Deng Fig-6B activation map.

Neural: Euler substeps at the template's native 0.1 ms inside the 2 ms
production cadence (20 substeps), physics at 1 ms (NSUB=2 per neural
step), 12 s runs. Perturbations: +2 nA at the neural port (depolarizing,
same sign as the twin) on RG_HC_flx or PF_ka_HC_flx, t = 2.0-2.1 s
(paper-exact) or 2.0-2.5 s (twin-window variant).

Outputs into this folder: deng_mujoco_<cond>.npz + a printed metrics
summary. Conditions: base / p1_rg / p2_pf / p3_rg_long.
"""
from pathlib import Path

import mujoco
import numpy as np

import sys

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))          # spinal/  -> deng_cpg_ode
import deng_cpg_ode as D                      # noqa: E402

XML = HERE.parent / "w2l_mujoco" / "w2l_air.xml"   # welded-root variant
DT_NET, DT_PHY = 0.002, 0.001
NSUB = int(DT_NET / DT_PHY)                   # physics steps per net step
NSUB_NEU = 20                                 # 0.1 ms Euler substeps
DUR = 12.0
CAP = 1.0

ACT = lambda v: np.clip(1.0 / (1.0 + np.exp(0.1532 * (-70.0 - v))) - 0.01,
                        0.0, 1.0)             # Deng Fig 6B map
MUS = {"hip_L_ext": D.NAMES.index("MN_hip_ext"),
       "hip_L_flx": D.NAMES.index("MN_hip_flx"),
       "knee_L_ext": D.NAMES.index("PF_ka_HC_ext"),
       "knee_L_flx": D.NAMES.index("PF_ka_HC_flx"),
       "ankle_L_ext": D.NAMES.index("PF_ka_HC_ext"),
       "ankle_L_flx": D.NAMES.index("PF_ka_HC_flx")}


def step_neural(V, h, iapp):
    """20 Euler substeps at 0.1 ms - verbatim inner loop of
    deng_cpg_ode.simulate (fixed-tau mode). NOTE: the ODE is
    milliseconds-native (nA/nF = mV/ms, tau_h in ms) - dt here is
    DT_NET in MILLISECONDS / NSUB_NEU (feeding seconds made every
    substep 1000x too small; caught by the probe script)."""
    dt = (DT_NET * 1000.0) / NSUB_NEU
    for _ in range(NSUB_NEU):
        sat = np.clip((V[D.src] - D.elo) / (D.ehi - D.elo), 0.0, 1.0)
        isyn = np.zeros(len(D.NAMES))
        np.add.at(isyn, D.dst, D.g * sat * (D.esyn - V[D.dst]))
        m_inf = 1 / (1 + D.KM * np.exp(D.SM * (D.EM - V)))
        h_inf = 1 / (1 + D.KH * np.exp(D.SH * (D.EH - V)))
        i_na = np.zeros(len(D.NAMES))
        na = list(D.NA)
        i_na[na] = D.G_NA * m_inf[na] * h[na] ** D.KH * (D.E_NA - V[na])
        V = V + (D.GM * (D.VREST - V) + iapp + isyn + i_na) * dt / D.CM
        h = h + dt * (h_inf - h) / D.TAU_H
    return V, h


def run(cond):
    pert_idx, t0w, t1w = {"base": (None, 0, 0),
                          "p1_rg": (D.NAMES.index("RG_HC_flx"), 2.0, 2.1),
                          "p2_pf": (D.NAMES.index("PF_ka_HC_flx"), 2.0, 2.1),
                          "p3_rg_long": (D.NAMES.index("RG_HC_flx"), 2.0, 2.5),
                          }[cond]
    m = mujoco.MjModel.from_xml_path(str(XML))
    d = mujoco.MjData(m)
    m.dof_damping[:] = 0.5 * m.dof_damping[:]   # gentler hang (documented
    # stand-in: full W2L damping is tuned for contact; we are in air)
    act_ids = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
               for a in MUS}
    V = D.VREST.copy()
    h = np.ones(len(D.NAMES)) * (1 / (1 + D.KH * np.exp(D.SH * (D.EH - D.VREST_HC))))
    h[[12, 13]] = 0.0
    nsteps = int(DUR / DT_NET)
    kick0, kick1 = 50, 51                      # t = 0.100-0.102 s (float-
    # safe integer indices; k*0.002 never == 0.1 exactly - the original
    # float-compare silently skipped the kick entirely)
    w0 = int(t0w / DT_NET) if pert_idx is not None else 0
    w1 = int(t1w / DT_NET) if pert_idx is not None else 0
    vE = np.zeros(nsteps + 1); vF = np.zeros(nsteps + 1)
    vE[0], vF[0] = V[0], V[1]
    qlog = np.zeros((nsteps + 1, 3))
    for k in range(nsteps):
        iapp = np.zeros(len(D.NAMES))
        if kick0 <= k < kick1:
            iapp[0] = 10.0                       # kickoff (template-verbatim)
        if pert_idx is not None and w0 <= k < w1:
            iapp[pert_idx] += 2.0                # Deng 2023 protocol
        V, h = step_neural(V, h, iapp)
        for a, ni in MUS.items():
            d.ctrl[act_ids[a]] = min(ACT(V[ni]), CAP)
        for _ in range(NSUB):
            mujoco.mj_step(m, d)
        vE[k + 1], vF[k + 1] = V[0], V[1]
        qadr_h = m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "hip_L")]
        qadr_k = m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "knee_L")]
        qadr_a = m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "ankle_L")]
        qlog[k + 1] = [d.qpos[qadr_h], d.qpos[qadr_k], d.qpos[qadr_a]]
    tt = np.arange(nsteps + 1) * DT_NET
    np.savez(HERE / f"deng_mujoco_{cond}.npz", t=tt, vE=vE, vF=vF, q=qlog)
    return tt, vE, vF


def onsets(tt, v, thr=-45.0):
    """Burst onsets: upward crossings of thr (mV) with 0.2 s refractory."""
    on = []
    last = -1e9
    above = v > thr
    for i in range(1, len(v)):
        if above[i] and not above[i - 1] and tt[i] - last > 0.2 and tt[i] > 1.0:
            on.append(tt[i]); last = tt[i]
    return np.array(on)


def report(tag, tt, vE, vF):
    oE, oF = onsets(tt, vE), onsets(tt, vF)
    pE_pre = np.mean(np.diff(oE[oE < 2.0])) if np.sum(oE < 2.0) >= 2 else np.nan
    pE_post = np.mean(np.diff(oE[oE > 2.6])) if np.sum(oE > 2.6) >= 2 else np.nan
    n_in = np.sum((oE >= 2.0) & (oE <= 2.6))
    print(f"[{tag}] ext onsets {len(oE)} (pre-per {pE_pre:.4f} s, in[2.0-2.6] "
          f"{n_in}, post-per {pE_post:.4f} s)")
    print(f"      next flx onsets around window: "
          f"{np.round(oF[(oF > 1.5) & (oF < 3.0)], 3)}")
    return pE_pre, pE_post, n_in, oF


if __name__ == "__main__":
    res = {}
    for cond in ("base", "p1_rg", "p2_pf", "p3_rg_long"):
        tt, vE, vF = run(cond)
        res[cond] = report(cond, tt, vE, vF)
    # phase-advance vs own baseline (wrapped to +- half period)
    base_f = res["base"][3]
    for cond in ("p1_rg", "p2_pf", "p3_rg_long"):
        pf = res[cond][3]
        if len(pf) and len(base_f):
            b_next = base_f[base_f >= 2.0][:1]
            p_next = pf[pf >= 2.0][:1]
            if len(b_next) and len(p_next):
                print(f"[{cond}] first flx onset at/after 2.0 s: "
                      f"{p_next[0]:.3f} vs baseline {b_next[0]:.3f} s "
                      f"(delta {p_next[0] - b_next[0]:+.3f} s)")
    print("DONE")
