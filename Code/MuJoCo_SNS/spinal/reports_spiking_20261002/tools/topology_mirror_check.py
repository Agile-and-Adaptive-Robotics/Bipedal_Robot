"""GATE 1 (SPIKING_MIRROR_PLAN.md): topology mirror check.

Builds the non-spiking network (build_network.py) and the spiking mirror
(build_network_spiking.py) under IDENTICAL gain configurations and
asserts the CONNECTION-CLASS CONTRACT is mirrored verbatim:
  - identical (source, destination, sign) edge multisets
  - identical input-port lists and MN/aff names
  - identical conditional-population presence (stage-1 gains: all
    conditional blocks ABSENT in both)
The ONLY allowed difference: the mirror's RD_* readout tap neurons
(pure sinks appended for the runner-facing analog levels) plus their
single tap edge each - counted and asserted exactly.

Configs probed (mirrors _fix_check.py):
  stage1 : all conditional gains 0
  tuned  : the _fix_check tuned set (every pathway nonzero)
  joint  : joint-layer PF flavor of the same tuned set

Usage: D:\\Anaconda\\envs\\myo\\python.exe topology_mirror_check.py
"""
import io
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

SPINAL = Path(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
os.chdir(SPINAL)
sys.path.insert(0, str(SPINAL))
assert "AARL_NET" not in os.environ

import mujoco

import params as P
import build_network as bn
import build_network_spiking as bs

ACTS = ["vas_lat_r", "semimem_r", "iliopsoas_r", "tib_ant_r",
        "vas_lat_l", "semimem_l", "iliopsoas_l", "tib_ant_l"]

STAGE1 = dict(phase_reset_e=0.0, phase_reset_f=0.0, f1_kneext_inh=0.0,
              f1_anklepf_inh=0.0, renshaw=0.0, ia_in=0.0, heel_rge=0.0,
              toe_rge=0.0, ib_rge=0.0)
TUNED = dict(phase_reset_e=0.5, phase_reset_f=0.5, f1_kneext_inh=0.6,
             f1_anklepf_inh=0.6, renshaw=0.5, ia_in=0.6, heel_rge=0.6,
             toe_rge=0.4, ib_rge=0.6, rg_weak_exc=0.0,
             ib_e_central=0.5, ia_f_central=0.5, ii_f_central=0.3,
             ii_e_central=0.3, ia_f_contra_f=0.4, v3_to_ibexc=0.4,
             c1_gain=0.6, v3_gain=0.25, full_rules=1.0,
             heel_pf_layer=0.5, toe_df_inh=1.0, heel_in_f_exc=0.5,
             ia_pf_f=0.5, ii_pf_f=0.5,
             aff_e_rg=0.3, aff_f_rg=0.3, aff_e_pf=0.3, aff_f_pf=0.3,
             vest_ext=0.3, vest_flex_inh=0.2, contra_kinh=0.5)


def edges_of(net):
    names = [p["name"] for p in net.net.populations]
    g = defaultdict(int)
    for c in net.net.connections:
        s, d = names[c["source"]], names[c["destination"]]
        sign = "exc" if c["params"].get("reversal_potential", 0) > -1e-6 \
            else "inh"
        g[(s, d, sign)] += 1
    return names, g


def run_config(label, gains, joint_pf):
    keep = dict(P.G)
    P.G.update(gains)
    P.G["joint_pf"] = 1.0 if joint_pf else 0.0
    m = mujoco.MjModel.from_xml_path(str(
        SPINAL.parents[2] / "Solid_Models" / "OpenSim" /
        "Gait2392_Robotbody" / "mjc" / "gait2392_simbody" /
        "gait2392_simbody_cvt3.xml"))
    acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(m.nu)]
    n_ns = bn.build(acts, interleg=True)
    n_sp = bs.build(acts, interleg=True)
    P.G.clear(); P.G.update(keep)

    names_ns, g_ns = edges_of(n_ns)
    names_sp, g_sp = edges_of(n_sp)

    # strip the mirror's readout taps (documented, counted)
    taps = [k for k in g_sp if k[1].startswith("RD_")]
    # per side: RG_E, RG_F + the PF cells that exist (4 phase / 6 joint)
    n_taps_expected = len(n_sp.sides) * (2 + (6 if joint_pf else 4))
    ok_taps = len(taps) == n_taps_expected
    for k in taps:
        del g_sp[k]
    rd_names = [n for n in names_sp if n.startswith("RD_")]
    ok_rd_names = len(rd_names) == n_taps_expected

    only_ns = {k: v for k, v in g_ns.items() if g_sp.get(k, 0) != v}
    only_sp = {k: v for k, v in g_sp.items() if g_ns.get(k, 0) != v}
    ok_edges = not only_ns and not only_sp
    ok_inputs = n_ns.inputs == n_sp.inputs
    ok_mn = n_ns.mn_names == n_sp.mn_names

    # stage-1 conditional absence (mirror of _fix_check (3))
    ok_cond = True
    if label == "stage1":
        bad = [n for n in names_ns if n.split("_")[0] in
               ("HEEL", "TOE", "LBIN", "IaIN", "RC", "PRESET", "PREA",
                "VEST", "AFF", "TOEDF", "IIX", "IBIN", "IIIN",
                "KINH")]
        ok_cond = not bad
        bad_sp = [n for n in names_sp if n.split("_")[0] in
                  ("HEEL", "TOE", "LBIN", "IaIN", "RC", "PRESET", "PREA",
                   "VEST", "AFF", "TOEDF", "IIX", "IBIN", "IIIN", "KINH")]
        ok_cond = ok_cond and not bad_sp

    print(f"[{label}] ns: {len(names_ns)} pops / {sum(g_ns.values())} "
          f"edges; sp: {len(names_sp)} pops / {sum(g_sp.values())} edges "
          f"(+{len(taps)} RD taps)")
    if only_ns:
        print(f"  ONLY-IN-NON-SPIKING: {list(only_ns.items())[:8]}")
    if only_sp:
        print(f"  ONLY-IN-SPIKING: {list(only_sp.items())[:8]}")
    print(f"  edges-mirror={ok_edges} inputs={ok_inputs} mn={ok_mn} "
          f"taps={ok_taps}({len(taps)}/{n_taps_expected}) "
          f"cond-absent={ok_cond}")
    ok = ok_edges and ok_inputs and ok_mn and ok_taps and ok_rd_names \
        and ok_cond
    print(f"  GATE1[{label}] {'PASS' if ok else 'FAIL'}")
    return ok


def main() -> int:
    ok = True
    ok &= run_config("stage1", STAGE1, joint_pf=False)
    ok &= run_config("tuned", TUNED, joint_pf=False)
    ok &= run_config("joint", TUNED, joint_pf=True)
    print("GATE1 TOPOLOGY MIRROR:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
