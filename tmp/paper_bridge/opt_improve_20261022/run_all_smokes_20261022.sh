#!/usr/bin/env bash
# Smoke-run driver for the 20261022 optimizer-improvement core deliverables.
# Sequential on purpose (10-core box; each MATLAB -batch is its own process).
set -u
MATLAB="D:/Program Files/MATLAB/R2025a/bin/matlab.exe"
OUT="D:/GitHub/Bipedal_Robot/tmp/paper_bridge/opt_improve_20261022"
REPO="D:/GitHub/Bipedal_Robot"
cd "$OUT"

echo "=== [1/5] BiPulley2 dry-run: deterministic replay gate over ALL 20 anchor mats ===" | tee -a run_all_smokes_20261022.log
"$MATLAB" -batch "run('Opt_run_BiPulley2_20261022.m')" > bi2_gate_dryrun_20261022.log 2>&1
echo "exit=$? (dry-run gate)" | tee -a run_all_smokes_20261022.log

echo "=== [2/5] BiPulley2 SMOKE mode=sweep (1 site, sub-grid {(1,1),(2,2)}, 60-eval budget) ===" | tee -a run_all_smokes_20261022.log
OPT_BPULLEY2_SMOKE=1 OPT_BPULLEY2_MODE=sweep "$MATLAB" -batch "run('Opt_run_BiPulley2_20261022.m')" > bi2_smoke_sweep_20261022.log 2>&1
echo "exit=$? (smoke sweep)" | tee -a run_all_smokes_20261022.log

echo "=== [3/5] BiPulley2 SMOKE mode=escalate (ladder seed (2,2), 2 rungs, tiny budget) ===" | tee -a run_all_smokes_20261022.log
OPT_BPULLEY2_SMOKE=1 OPT_BPULLEY2_MODE=escalate "$MATLAB" -batch "run('Opt_run_BiPulley2_20261022.m')" > bi2_smoke_escalate_20261022.log 2>&1
echo "exit=$? (smoke escalate)" | tee -a run_all_smokes_20261022.log

echo "=== [4/5] Opt_sanity_BiPulley (Ben's sanity gate, read-only run: G=1 identity 1e-8 / G=2 closure 1e-6) ===" | tee -a run_all_smokes_20261022.log
cd "$REPO/Code/Matlab/Mesh_Optimization"
"$MATLAB" -batch "run('Opt_sanity_BiPulley.m')" > "$OUT/bi_sanity_BiPulley_20261022.log" 2>&1
echo "exit=$? (sanity BiPulley)" | tee -a "$OUT/run_all_smokes_20261022.log"

echo "=== [5/5] Flexor sweep33 SMOKE (OPT_PULLEY_SMOKE=1, sub-grid nB=1 G=2, embedded sanity checks) ===" | tee -a "$OUT/run_all_smokes_20261022.log"
cd "$OUT"
OPT_PULLEY_SMOKE=1 "$MATLAB" -batch "run('Opt_run_pulley_sweep33_20261022.m')" > flexor_sweep33_smoke_20261022.log 2>&1
echo "exit=$? (flexor sweep33 smoke)" | tee -a run_all_smokes_20261022.log

echo "ALL SMOKES DONE $(date)" | tee -a run_all_smokes_20261022.log
