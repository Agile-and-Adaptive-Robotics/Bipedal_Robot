r"""Write w2l_walk_best_params.json from the winner probe npz (egg,
campaign 2026-10-16, stage 'bootstrap driver (robot)')."""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DT = 0.001
z = np.load(os.path.join(HERE, "probe_w01.npz"))
t = z["t"]
w = t >= 4.5
tw = t[w]
heel = z["heelF"][w]
rg = z["rgE"][w]
com = z["com"][w]


def primary(f, thr=20.0, merge_s=0.45):
    on = f > thr
    st = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    prim = []
    for x in st:
        if not prim or (tw[x] - tw[prim[-1]]) > merge_s:
            prim.append(int(x))
        # double-tap follower: same cycle
    return prim, [int(x) for x in st]


pL, sL = primary(heel[:, 0])
pR, sR = primary(heel[:, 1])
perL = float(np.median(np.diff(tw[pL]))) if len(pL) >= 3 else float("nan")
perR = float(np.median(np.diff(tw[pR]))) if len(pR) >= 3 else float("nan")
# inter-leg phase from primary strikes
ph = []
for a in range(len(pR) - 1):
    for b in pL:
        if tw[pR[a]] <= tw[b] <= tw[pR[a + 1]]:
            ph.append((tw[b] - tw[pR[a]]) /
                      (tw[pR[a + 1]] - tw[pR[a]]))
r = float(np.corrcoef(rg[::5, 0], rg[::5, 1])[0, 1])
fwd = float(com[-1, 0] - com[0, 0])
duty_h = [float((heel[:, k] > 5).mean()) for k in (0, 1)]
toe = z["toeF"][w]
duty_c = [float(((heel[:, k] > 5) | (toe[:, k] > 5)).mean()) for k in (0, 1)]
# sustained-phase speed: linear fit of com x over t in [8, 30] s
m2 = tw >= 8.0
slope = float(np.polyfit(tw[m2], com[m2, 0], 1)[0])

out = {
  "winner": {
    "model": "W2L robot, w2l_mujoco (w2l_mjcf_fixed.xml body + M5 "
             "afferented split-RG net, build_w2l_aff_net.build)",
    "knobs": {
      "te": 4.0, "tf": 5.0,
      "te_l": None, "te_r": None, "tf_l": None, "tf_r": None,
      "alpha_l": 0.0, "alpha_r": 0.0,
      "legcoll": 0, "sink": 0.06, "ctn_amp": 2.0,
      "rig": 0.25, "walker_yield": 1.5, "dur_walk": 30.0
    },
    "protocol_defaults": {"hold": 3.0, "lower": 1.0, "wsup": 0.7,
                          "drop": 0.026, "cap": 0.5, "acap": -1.0,
                          "cref": 50.0, "ibnA": 1.0, "tau": 0.25,
                          "damp": 3.0, "settle": 0.4, "comm": 1.0},
    "command": ("cd /d D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal"
                "\\w2l_mujoco && D:\\Anaconda\\envs\\myo\\python.exe "
                "test_w2l_ground.py --phase=walk --te=4 --tf=5 "
                "--legcoll=0 --sink=0.06 --ctn_amp=2 --rig=0.25 "
                "--walker-yield=1.5 --dur-walk=30"),
    "metrics_25p5s_window": {
      "heel_strikes_L": len(sL), "heel_strikes_R": len(sR),
      "primary_stride_cycles_L": len(pL), "primary_stride_cycles_R": len(pR),
      "period_L_s": round(perL, 3), "period_R_s": round(perR, 3),
      "rg_e_corr": round(r, 3),
      "interleg_phase_median": round(float(np.median(ph)), 3),
      "forward_disp_m": round(fwd, 3),
      "mean_forward_speed_mps": round(fwd / float(tw[-1] - tw[0]), 4),
      "sustained_speed_fit_8_30s_mps": round(slope, 4),
      "heel_duty_L": round(duty_h[0], 3), "heel_duty_R": round(duty_h[1], 3),
      "combined_duty_L": round(duty_c[0], 3),
      "combined_duty_R": round(duty_c[1], 3),
      "fall": False, "finite": True,
      "tilt_max_deg": 16.1, "pelvis_z_min_m": 0.85
    }
  },
  "runner_up_more_upright": {
    "knobs": {"te": 4.0, "tf": 5.0, "legcoll": 0, "sink": 0.06,
              "ctn_amp": 2.0, "rig": 0.5, "walker_yield": 3.0,
              "dur_walk": 30.0},
    "metrics": {"strikes": "24/37", "rg_e_corr": -0.770,
                "fwd_m_25p5s": 0.437, "tilt_max_deg": 7.3,
                "heel_duty": "0.065/0.068"}
  },
  "provenance": {
    "stage": "bootstrap driver (robot), overnight campaign 2026-10-16",
    "sweep_runs": 128,
    "sweeps": "A symmetric te x tf (25) | B E-dominant+asym gamma (20) | "
              "C alpha (12) | D stance/harness (14) | E legcoll+sink (13) "
              "| F refine (18) | G rig/progression (10) | H rig refine (8) "
              "| W yielding walker (8)",
    "results_csv": "reports_20261016/w2l_walk_sweep/results.csv",
    "probe_npz": "reports_20261016/w2l_walk_sweep/probe_w01.npz",
    "report": "reports_20261016/W2L_ROBOT_BOOTSWEEP_20261016.md"
  }
}
json.dump(out, open(os.path.join(HERE, "w2l_walk_best_params.json"), "w"),
          indent=1)
print("primary L cycles:", len(pL), "R:", len(pR))
print("periods:", round(perL, 3), round(perR, 3))
print("interleg phase median:", round(float(np.median(ph)), 3),
      "n =", len(ph))
print("sustained slope m/s:", round(slope, 4))
print("wrote w2l_walk_best_params.json")
