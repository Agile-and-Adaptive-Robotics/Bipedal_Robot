#!/bin/bash
# BiPulley expansion campaign launch, easteregg2, 2026-09-26.
# Sanity gates PASSED earlier this session (sanity_pulley_20260926.log,
# sanity_BiPulley_20260926.log in this folder). This script launches the two
# production runs SEQUENTIALLY (one parpool(6) at a time on the 10-core box):
#   1. Opt_run_pulley   -- cwd = Testing_Data/2022_02_Festo (cwd gotcha:
#      buildKneeFlexorContext20mm.m:129 readmatrix('OpenSim_Bifem_Results.txt')
#      resolves against the CURRENT FOLDER; the file exists only there).
#      Launched BY NAME after addpath (runs in BASE; no run('fullpath')).
#      In-file settings kept verbatim: OPT_PULLEY_SMOKE unset (full run,
#      liveRun auto-armed -> dated full-workspace save), OPT_PULLEY_CONFIGS
#      unset (single config nPulleyBPA=2, G=1, moving_via), surrogateopt
#      1000 + patternsearch 15000 evals, parpool(6) opened in-script at
#      Opt_run_pulley.m:174.
#   2. Opt_run_BiPulley -- cwd = Mesh_Optimization; RUN_BATCH=true +
#      FULL_RUN=true (flipped in-file 2026-09-26, budgets untouched:
#      surrogateopt 1000 + patternsearch 5000 x 5 OpenSim muscles);
#      parpool(6) PRE-OPENED here so the bare `parpool` at
#      Opt_run_BiPulley.m:420 is skipped by its isempty(gcp('nocreate'))
#      guard (cap 6 per campaign directive; default would be 10).
unset OPT_PULLEY_SMOKE OPT_PULLEY_CONFIGS OPT_BIPULLEY_SMOKE
DIG="D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo/Dig_out/pulley_campaign_20260926"
MATLAB="D:/Program Files/MATLAB/R2025a/bin/matlab.exe"
MO="D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization"
FESTO="D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo"

LAUNCH1="matlab -batch \"addpath('$MO'); Opt_run_pulley\"   (cwd = $FESTO)"
LAUNCH2="matlab -batch \"parpool(6); addpath('$MO'); Opt_run_BiPulley\"   (cwd = $MO)"
{
  echo "Campaign launch record 2026-09-26"
  echo "hostname: $(hostname)"
  echo "env guard: OPT_PULLEY_SMOKE / OPT_PULLEY_CONFIGS / OPT_BIPULLEY_SMOKE unset"
  echo "LAUNCH1 (Opt_run_pulley):    $LAUNCH1"
  echo "LAUNCH2 (Opt_run_BiPulley):  $LAUNCH2"
} > "$DIG/launch_record_20260926.txt"

echo "[chain] START Opt_run_pulley $(date +%Y-%m-%dT%H:%M:%S) :: $LAUNCH1"
cd "$FESTO" || exit 1
"$MATLAB" -batch "addpath('$MO'); Opt_run_pulley" > "$DIG/prod_Opt_run_pulley_20260926.log" 2>&1
echo "[chain] EXIT_PULLEY=$? $(date +%Y-%m-%dT%H:%M:%S)"

echo "[chain] START Opt_run_BiPulley $(date +%Y-%m-%dT%H:%M:%S) :: $LAUNCH2"
cd "$MO" || exit 1
"$MATLAB" -batch "parpool(6); addpath('$MO'); Opt_run_BiPulley" > "$DIG/prod_Opt_run_BiPulley_20260926.log" 2>&1
echo "[chain] EXIT_BIPULLEY=$? $(date +%Y-%m-%dT%H:%M:%S)"
echo "[chain] DONE $(date +%Y-%m-%dT%H:%M:%S)"
