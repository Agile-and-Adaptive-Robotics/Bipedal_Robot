"""Reflex-loop sign audit (2026-10-21, Ben's REFLEXIVE-balance directive).

Dynamic probe of the Ia/II/Ib arcs at the QUIET-STANDING operating point
(DRIVE=0 nA, POSTURE=1 nA, the drive_posture schedule during stand1/stand2),
per model family. For each arc: inject a unit port current into one
afferent encoder, step to steady state, read the MN response.

Verdict per arc: does Ia EXCITE the homonymous MN (monosynaptic,
anti-stretch negative feedback) and INHIBIT the antagonist; does II excite
homonymous; does Ib autogenically inhibit (and is the stance IBEXC reversal
path alive at DRIVE=0, i.e. is its RG_E gate open?).

Also reads the runner's actual encoder normalizations at the standing pose
(len_norm per muscle from MuJoCo) to see whether the II length term has any
signal to work with at posture.

Families probed: stock (s3k default build) and s3k-trial34 topology
(full_rules=1, ia_in=0.825) - the two the anchor/winner campaigns use.
"""
from __future__ import annotations

import numpy as np
import mujoco

import build_network as bn
import params as P
from params import DT, E_HI, G
import runner as R


def steady(net, u, seconds=1.5):
    n = int(seconds / DT)
    for _ in range(n):
        net.step(u)
    return net.compiled.V.copy()


def probe_family(tag, acts, gain_overrides=None):
    # mutate THE dict build_network bound at import (reload(params) would
    # orphan it - the first run's second family silently probed stock twice)
    saved = {k: G[k] for k in (gain_overrides or {})}
    for k, v in (gain_overrides or {}).items():
        G[k] = v
    net = bn.build(acts, dt=DT)
    for k, v in saved.items():
        G[k] = v
    iport = {p: i for i, p in enumerate(net.inputs)}

    # --- quiet-standing operating point: DRIVE 0, POSTURE 1, no afferents
    u = net.make_inputs()
    u[iport["DRIVE"]] = 0.0
    u[iport["POSTURE"]] = 1.0
    v0 = steady(net, u)
    rge = v0[net.idx["RG_E_r"]]
    rgf = v0[net.idx["RG_F_r"]]
    stance_gate = float(np.clip(rge / E_HI, 0, 1))
    print(f"\n== {tag} ==")
    print(f"steady standing (DRIVE=0,POSTURE=1): RG_E_r={rge:.3f} mV "
          f"RG_F_r={rgf:.3f} mV -> stance gate={stance_gate:.3f}")
    print(f"  -> effective extensor gains at this gate: "
          f"Ia x{0.5 + 0.5 * stance_gate:.2f}, II x{0.3 + 0.7 * stance_gate:.2f}, "
          f"Ib x{1.0 - 0.7 * stance_gate:.2f} "
          f"(IBEXC gate conductance 1.0 x RG_E act {stance_gate:.2f})")

    def arc(port, current, mns, label):
        u2 = u.copy()
        u2[iport[port]] = current
        net.step(u2)
        v1 = steady(net, u2, seconds=0.5)
        net.step(u)  # settle back partially
        base = steady(net, u, seconds=0.5)
        out = {m: float(v1[net.idx[f"MN_{m}"]] - base[net.idx[f"MN_{m}"]])
               for m in mns}
        print(f"  {label:46s} " + "  ".join(
            f"dMN_{m}={d:+.3f} mV" for m, d in out.items()))
        return out

    # --- soleus_r (ankle_pf, stance extensor) arcs
    arc("Ia_soleus_r", 1.6, ["soleus_r", "tib_ant_r"],
        "Ia_soleus_r +1.6nA (stretch):")
    arc("II_soleus_r", 1.2, ["soleus_r", "tib_ant_r"],
        "II_soleus_r +1.2nA (length):")
    arc("Ib_soleus_r", 1.5, ["soleus_r", "tib_ant_r"],
        "Ib_soleus_r +1.5nA (force):")
    # --- tib_ant_r (ankle_df, antagonist) arcs
    arc("Ia_tib_ant_r", 1.6, ["tib_ant_r", "soleus_r"],
        "Ia_tib_ant_r +1.6nA (stretch):")
    # --- vastus (knee_ext stance group): the IBEXC group member
    arc("Ib_vas_lat_r", 1.5, ["vas_lat_r", "bifemsh_r"],
        "Ib_vas_lat_r +1.5nA (force):")

    # --- Ib net effect vs stance: repeat soleus Ib with DRIVE raised
    u3 = u.copy()
    u3[iport["DRIVE"]] = 4.0
    v_hi = steady(net, u3)
    rge_hi = v_hi[net.idx["RG_E_r"]]
    gate_hi = float(np.clip(rge_hi / E_HI, 0, 1))
    u3[iport[f"Ib_soleus_r"]] = 1.5
    v1 = steady(net, u3, seconds=0.5)
    d_hi = float(v1[net.idx["MN_soleus_r"]] - v_hi[net.idx["MN_soleus_r"]])
    ibexc_v = float(v1[net.idx.get(f"IBEXC_ankle_pf_r", net.idx["RG_E_r"])]) \
        if f"IBEXC_ankle_pf_r" in net.idx else float("nan")
    print(f"  DRIVE=4: RG_E_r={rge_hi:.2f} (gate {gate_hi:.2f}); "
          f"Ib_soleus_r +1.5nA -> dMN_soleus={d_hi:+.3f} mV "
          f"(load-sharing reversal {'ALIVE' if d_hi > 0 else 'DEAD/negative'}; "
          f"IBEXC_ankle_pf_r V={ibexc_v:.2f} mV)")
    return net


def standing_pose_signals():
    """Runner encoder normalizations at the standing pose (len_norm sign
    for the antigravity muscles: is there II length signal at posture?)."""
    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    R.apply_start_pose(model, data)
    key = R.capture_pose(model, data)
    model = R.apply_harness(model, data, kxy=1500.0)
    data = mujoco.MjData(model)
    R.seed_pose(model, data, key)
    mujoco.mj_forward(model, data)
    Lrange = model.actuator_lengthrange
    Lmid, Lhalf = Lrange.mean(axis=1), np.maximum(
        (Lrange[:, 1] - Lrange[:, 0]) / 2, 1e-3)
    L = data.actuator_length
    len_norm = np.clip((L - Lmid) / Lhalf, -1.0, 1.0)
    names = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
             for i in range(model.nu)]
    print("\n== standing-pose II length signals (len_norm = (L-Lmid)/Lhalf) ==")
    for m in ("soleus_r", "med_gas_r", "vas_lat_r", "vas_med_r", "glut_max2_r",
              "ercspn_r", "tib_ant_r", "rect_fem_r", "semimem_r", "iliacus_r"):
        i = names.index(m)
        print(f"  {m:14s} L={L[i]:.4f} mid={Lmid[i]:.4f} "
              f"range=[{Lrange[i,0]:.4f},{Lrange[i,1]:.4f}] "
              f"len_norm={len_norm[i]:+.3f}")


def main():
    acts = None
    # actuator list from the model
    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(model.nu)]
    probe_family("STOCK default build (s3k anchor flavor, full_rules=0)", acts)
    probe_family(
        "S3K trial-34 topology (full_rules=1, ia_in=0.825, ib gains on)",
        acts,
        gain_overrides=dict(full_rules=1.0, ia_in=0.8250460785281135,
                            ib_e_central=0.24354654619098992,
                            ib_rge=0.6741474252269108,
                            heel_rge=0.37530421840765804,
                            toe_rge=0.14572281192093522))
    standing_pose_signals()


if __name__ == "__main__":
    main()
