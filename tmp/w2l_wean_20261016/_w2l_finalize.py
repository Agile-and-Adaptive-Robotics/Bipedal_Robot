r"""FINALIZE the W2L wean-stage winner param file (2026-10-16).

Authoritative assembly after the CLI-reproduction audit: the in-process
campaign verify at rungs 0.50-0.20 turned out context-dependent (the
r020 'survivor' + 5/12 basin failures were the tell); the ROBUST frontier
is the fresh-process CLI reproduction of every rung. Winner = the deepest
ROBUST rung (0.60). Also patches ladder_summary.json with the CLI view.
"""
from __future__ import annotations

import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

BOOT_TE, BOOT_TF = 4.0, 5.0
BASE = dict(legcoll=0, sink=0.06, rig=0.25, walker_yield=1.5)


def load(name):
    return json.load(open(os.path.join(HERE, name), encoding="utf-8"))


def main() -> int:
    lad = load("ladder_summary.json")
    cli = load("cli_verify_rungs.json")
    robust = [c["dscale"] for c in cli if c["cli_survive"]]
    dead = [c["dscale"] for c in cli if not c["cli_survive"]]
    winner_ds = min(robust)
    win = load(f"round_best_r{int(round(winner_ds*100)):03d}.json")
    wcli = next(c for c in cli if abs(c["dscale"] - winner_ds) < 1e-9)
    p = win["params"]
    bas60 = load("basin_results_r060.json")
    bas20 = load("basin_results_r020.json")
    bas100 = load("basin_results_r100.json")

    # baseline row (bootstrap config, CLI context) from the trials CSV
    base_row = None
    with open(os.path.join(HERE, "wean_trials.csv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row.get("tag") == "baseline_boot_d30":
                base_row = row
    assert base_row is not None, "baseline_boot_d30 row missing"

    ds = winner_ds
    te, tf = BOOT_TE * ds, BOOT_TF * ds
    m = wcli["cli_metrics"]
    knobs = dict(
        te=round(te, 4), tf=round(tf, 4), dscale=ds,
        g_heel=p["g_heel"], g_ib=p["g_ib"], ibnA=p["ibnA"],
        ctn_amp=p["ctn_amp"], cref=p["cref"], comm=p["comm"],
        g_toe_in=round(1.0 * p["toe_m"], 4),
        g_toe_df=round(5.0 * p["toe_m"], 4),
        g_toe_df_mn=round(2.749 * p["toe_m"], 4),
        **BASE, dur_walk=30.0)

    ev_cmd = ("cd /d D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal"
              "\\reports_20261016\\w2l_wean && D:\\Anaconda\\envs\\myo"
              "\\python.exe _w2l_wean_eval.py --id=reproduce_finalist "
              f"--dur=30 --dscale={ds} --g-heel={p['g_heel']} "
              f"--g-ib={p['g_ib']} --ibnA={p['ibnA']} "
              f"--ctn-amp={p['ctn_amp']} --cref={p['cref']} "
              f"--comm={p['comm']} --toe-m={p['toe_m']}")
    cli_cmd = ("cd /d D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal"
               "\\w2l_mujoco && D:\\Anaconda\\envs\\myo\\python.exe "
               "test_w2l_ground.py --phase=walk "
               f"--te={te!r} --tf={tf!r} "
               f"--legcoll=0 --sink=0.06 --rig=0.25 --walker-yield=1.5 "
               f"--ctn-amp={p['ctn_amp']!r} --cref={p['cref']!r} "
               f"--ibnA={p['ibnA']!r} --dur-walk=30 "
               f"--g-heel={p['g_heel']!r} --g-ib={p['g_ib']!r} "
               f"--comm={p['comm']!r} "
               f"--g-toe-in={1.0*p['toe_m']!r} "
               f"--g-toe-df={5.0*p['toe_m']!r} "
               f"--g-toe-df-mn={2.749*p['toe_m']!r}")

    out = dict(
        winner=dict(
            model=("W2L robot, w2l_mujoco (w2l_mjcf_fixed.xml body + M5 "
                   "afferented split-RG net via build_w2l_aff_net), "
                   "drop-protocol harness + yielding virtual walker "
                   "(wsup 0.7, rig 0.25, walker_yield 1.5 s)"),
            stage=("feedback optimizer (robot): optuna feedback tuning "
                   "(study curr_s6_drive_robot_r060) + drive-wean ladder "
                   "+ CLI-reproduction audit + basin gate"),
            knobs=knobs,
            metrics_30s_cli_reproduced=m,
            metrics_30s_campaign_context=win["m30"],
            reproduction_note=("CLI fresh-process runs are deterministic "
                               "among themselves (4x identical, obj "
                               "94.86191); the campaign in-process verify "
                               "gives obj 100.64 -- cross-context BLAS "
                               "summation-order chaos, documented stack "
                               "property; ALL conclusions use the CLI "
                               "context"),
            basin="PASS 12/12 (+-1% joint perturbations, min cyc 20, "
                  "min fwd +2.145 m)",
            reproduce_evaluator_cmd=ev_cmd,
            reproduce_test_cli_cmd=cli_cmd),
        baseline_bootstrap=dict(
            note=("drive-only bootstrap winner (Ben's rules gains, te=4/"
                  "tf=5) same CLI context, 30 s"),
            metrics=dict(cyc=f"{base_row['cycL']}/{base_row['cycR']}",
                         phase_med=float(base_row["phase_med"]),
                         slope=float(base_row["slope"]),
                         duty_combined=(float(base_row["dutycL"]),
                                        float(base_row["dutycR"])),
                         tilt_max=float(base_row["tilt_max"]),
                         obj=float(base_row["obj"]))),
        wean_ladder_cli_robust=dict(
            robust_rungs=robust,
            dead_rungs=dead,
            death_scale=max(dead),
            winner_dscale=winner_ds,
            failure_mode_at_050=("inter-leg phase 0.322 (out of 0.35-0.65 "
                                 "band), cycR 9 -- antiphase breaks first"),
            failure_mode_below_040=("0 heel strikes either side, R foot "
                                    "planted (combined duty 1.0) -- gait "
                                    "dead"),
            marginal_in_process_datum=dict(
                dscale=0.20,
                what_happened=("in-process campaign verify PASSED the "
                               "survival criterion (34/33 cycles, 0.065 "
                               "m/s) but the SAME config in a fresh CLI "
                               "process is DEAD (0/0 cycles) and FAILED "
                               "the +-1% basin gate (5/12 die) -- a "
                               "context-dependent marginal basin, kept "
                               "as the cautionary datum"),
                basin_pass=False,
                sensitivity={k: v["obj"]
                             for k, v in bas20["sensitivity"].items()}),
            ladder_table=[dict(
                dscale=c["dscale"],
                cli_survive=c["cli_survive"],
                cyc=f"{c['cli_metrics'].get('cycL')}/"
                    f"{c['cli_metrics'].get('cycR')}",
                phase=c["cli_metrics"].get("phase_med"),
                slope=c["cli_metrics"].get("slope"),
                obj=c["cli_metrics"].get("obj")) for c in cli]),
        basin_and_sensitivity=dict(
            r060_basin=dict(pass_=bas60["basin_pass"],
                            min_cyc=bas60["basin_margin_cyc"],
                            min_fwd=bas60["basin_margin_fwd"],
                            sensitivity={k: dict(
                                obj=v["obj"], survive=v["survive"])
                                for k, v in bas60["sensitivity"].items()}),
            r100_basin=dict(pass_=bas100["basin_pass"],
                            min_cyc=bas100["basin_margin_cyc"],
                            min_fwd=bas100["basin_margin_fwd"],
                            sensitivity={k: dict(
                                obj=v["obj"], survive=v["survive"])
                                for k, v in bas100["sensitivity"].items()}),
            interpretation=(
                "r060: no single pathway is solely load-bearing (any one "
                "removed: obj 86-101 vs anchor 101) but Ben's UNTUNED "
                "default gains at 0.6 drive are DEAD (obj 7, cyc 10/0) -- "
                "the tuned SET is what sustains the weaned gait. r100 "
                "(full drive): defaults survive (94, the bootstrap), comm "
                "silencing hurts most (40). r020 (marginal): every "
                "pathway critical (5-7).")),
        provenance=dict(
            stage="feedback optimizer (robot), overnight campaign 2026-10-16",
            db="reports_20261016/w2l_wean/optuna_w2l_robot.db",
            studies=[f"curr_s6_drive_robot_r{int(round(r['dscale']*100)):03d}"
                     for r in lad["rounds"]],
            trials_csv="reports_2026106/w2l_wean/wean_trials.csv",
            cli_verify="reports_20261016/w2l_wean/cli_verify_rungs.json",
            basins=["basin_results_r100.json", "basin_results_r060.json",
                    "basin_results_r020.json"],
            bootstrap_report="reports_20261016/W2L_ROBOT_BOOTSWEEP_20261016.md"),
    )
    out["provenance"]["trials_csv"] = \
        "reports_20261016/w2l_wean/wean_trials.csv"
    path = os.path.join(HERE, "w2l_wean_best_params.json")
    json.dump(out, open(path, "w", encoding="utf-8"), indent=1)
    print("wrote", path)
    print("WINNER dscale", winner_ds, "CLI:", {k: m.get(k) for k in
          ("cycL", "cycR", "phase_med", "slope", "dutycL", "dutycR",
           "tilt_max", "obj", "survive")})

    # patch ladder_summary with the CLI-robust authoritative view
    lad["cli_robust_view"] = dict(robust_rungs=robust, dead_rungs=dead,
                                  death_scale=max(dead),
                                  winner_dscale=winner_ds)
    json.dump(lad, open(os.path.join(HERE, "ladder_summary.json"), "w",
                        encoding="utf-8"), indent=1)
    print("patched ladder_summary.json (cli_robust_view)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
