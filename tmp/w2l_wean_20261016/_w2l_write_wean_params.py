r"""Write w2l_wean_best_params.json — the WEAN-stage winner param file
(overnight campaign 2026-10-16, egg; stage "feedback optimizer (robot)").

Assembles: ladder_summary.json (finalist round) + basin_results.json
(basin gate + sensitivity) -> one canonical json with the FULL knob set,
the 30 s verification metrics, both reproduce commands, provenance.

Usage:
    D:\Anaconda\envs\myo\python.exe _w2l_write_wean_params.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

BOOT_TE, BOOT_TF = 4.0, 5.0
BASE = dict(legcoll=0, sink=0.06, rig=0.25, walker_yield=1.5)


def main() -> int:
    fin_name = sys.argv[1] if len(sys.argv) > 1 else "round_best_r060.json"
    bas_name = sys.argv[2] if len(sys.argv) > 2 else "basin_results_r060.json"
    lad = json.load(open(os.path.join(HERE, "ladder_summary.json"),
                         encoding="utf-8"))
    fin = json.load(open(os.path.join(HERE, fin_name), encoding="utf-8"))
    fin_ds = float(fin["dscale"])
    p = fin["params"]
    m = fin["m30"]
    bas = json.load(open(os.path.join(HERE, bas_name), encoding="utf-8"))
    # the deepest-wean datum (r020): survived the ask's criterion but
    # FAILED the +-1% basin gate -- recorded, flagged, not the winner
    marg = dict(dscale=0.20)
    try:
        m20 = json.load(open(os.path.join(HERE, "round_best_r020.json"),
                             encoding="utf-8"))
        b20 = json.load(open(os.path.join(HERE, "basin_results_r020.json"),
                             encoding="utf-8"))
        marg.update(
            params=m20["params"],
            m30_digest=dict(cyc=f"{m20['m30']['cycL']}/{m20['m30']['cycR']}",
                            phase=m20["m30"]["phase_med"],
                            slope=m20["m30"]["slope"],
                            fwd=m20["m30"]["fwd"],
                            tilt=m20["m30"]["tilt_max"],
                            obj=m20["m30"]["obj"],
                            survive=m20["m30"]["survive"]),
            basin_pass=b20["basin_pass"],
            basin_note=("5/12 joint +-1% perturbations die -- marginal "
                        "basin, flagged per campaign contracts"),
            sensitivity={k: dict(obj=v["obj"], survive=v["survive"])
                         for k, v in b20["sensitivity"].items()})
    except FileNotFoundError:
        marg["note"] = "r020 files not found"
    ds = fin_ds
    te, tf = BOOT_TE * ds, BOOT_TF * ds
    toe_m = p["toe_m"]
    knobs = dict(
        te=round(te, 4), tf=round(tf, 4), dscale=ds,
        g_heel=p["g_heel"], g_ib=p["g_ib"], ibnA=p["ibnA"],
        ctn_amp=p["ctn_amp"], cref=p["cref"], comm=p["comm"],
        g_toe_in=round(1.0 * toe_m, 4), g_toe_df=round(5.0 * toe_m, 4),
        g_toe_df_mn=round(2.749 * toe_m, 4),
        **BASE, dur_walk=30.0,
    )
    ev_cmd = ("cd /d D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal"
              "\\reports_20261016\\w2l_wean && D:\\Anaconda\\envs\\myo"
              "\\python.exe _w2l_wean_eval.py --id=reproduce_finalist "
              f"--dur=30 --dscale={ds} "
              f"--g-heel={p['g_heel']} --g-ib={p['g_ib']} "
              f"--ibnA={p['ibnA']} --ctn-amp={p['ctn_amp']} "
              f"--cref={p['cref']} --comm={p['comm']} "
              f"--toe-m={toe_m}")
    cli_cmd = ("cd /d D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal"
               "\\w2l_mujoco && D:\\Anaconda\\envs\\myo\\python.exe "
               "test_w2l_ground.py --phase=walk "
               f"--te={te:.4g} --tf={tf:.4g} "
               f"--legcoll=0 --sink=0.06 --rig=0.25 --walker-yield=1.5 "
               f"--ctn-amp={p['ctn_amp']} --cref={p['cref']} "
               f"--ibnA={p['ibnA']} --dur-walk=30 "
               f"--g-heel={p['g_heel']} --g-ib={p['g_ib']} "
               f"--comm={p['comm']} "
               f"--g-toe-in={1.0*toe_m:.4g} --g-toe-df={5.0*toe_m:.4g} "
               f"--g-toe-df-mn={2.749*toe_m:.4g}")
    out = dict(
        winner=dict(
            model=("W2L robot, w2l_mujoco (w2l_mjcf_fixed.xml body + M5 "
                   "afferented split-RG net via build_w2l_aff_net), "
                   "drop-protocol harness + yielding virtual walker"),
            stage=("feedback optimizer (robot): optuna feedback tuning + "
                   "drive-wean ladder, study curr_s6_drive_robot_r*"),
            knobs=knobs,
            metrics_30s=m,
            reproduce_evaluator_cmd=ev_cmd,
            reproduce_test_cli_cmd=cli_cmd),
        ladder=dict(
            rounds=[dict(dscale=r["dscale"], obj20=r["obj20"],
                         survive30=r["survive30"],
                         m30_digest=dict(
                             cyc=f"{r['m30']['cycL']}/{r['m30']['cycR']}",
                             phase=r["m30"]["phase_med"],
                             slope=r["m30"]["slope"],
                             fwd=r["m30"]["fwd"],
                             tilt=r["m30"]["tilt_max"],
                             obj=r["m30"]["obj"],
                             survive=r["m30"]["survive"]))
                    for r in lad["rounds"]],
            finalist_dscale=fin_ds,
            death_scale=lad.get("death_scale"),
            wall_min=lad.get("wall_min")),
        basin_and_sensitivity=dict(
            basin_pass=bas["basin_pass"],
            basin_margin_cyc=bas["basin_margin_cyc"],
            basin_margin_fwd=bas["basin_margin_fwd"],
            anchor_obj=bas["anchor"]["obj"],
            sensitivity={k: dict(obj=v["obj"], cyc=f"{v['cycL']}/{v['cycR']}",
                                 slope=v["slope"], fwd=v["fwd"],
                                 survive=v["survive"])
                         for k, v in bas["sensitivity"].items()}),
        marginal_deepest_wean=marg,
        provenance=dict(
            stage="feedback optimizer (robot), overnight campaign 2026-10-16",
            db="reports_20261016/w2l_wean/optuna_w2l_robot.db",
            trials_csv="reports_20261016/w2l_wean/wean_trials.csv",
            ladder="reports_20261016/w2l_wean/ladder_summary.json",
            basin="reports_20261016/w2l_wean/basin_results.json",
            bootstrap_report="reports_20261016/W2L_ROBOT_BOOTSWEEP_20261016.md"),
    )
    path = os.path.join(HERE, "w2l_wean_best_params.json")
    json.dump(out, open(path, "w", encoding="utf-8"), indent=1)
    print("wrote", path)
    print("finalist dscale", fin_ds, "| m30:", {k: m[k] for k in
          ("cycL", "cycR", "phase_med", "slope", "fwd", "tilt_max",
           "survive", "obj")})
    return 0


if __name__ == "__main__":
    sys.exit(main())
