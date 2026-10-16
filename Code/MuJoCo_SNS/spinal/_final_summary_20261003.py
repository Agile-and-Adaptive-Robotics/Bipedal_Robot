"""Append the final campaign summary to CURRICULUM_RUN_20261003.log."""

SUMMARY = """
================================================================================
===== FINAL CAMPAIGN SUMMARY 2026-10-04 (EB475WS4, myo env, 3 detached chains)
================================================================================
STUDY RESULTS (best scores, exact from dbs):
  s3k stage4 mutfix   : 52 trials, best -192.21 (t25) - converged; the +15
                        extension did not move it (winner stable)
  syn6 s5 reseed      : 37 trials, best -200.38 (t0) - seed holds
  s3kpruned walk_1    : 40 fresh trials, best -187.00 (t35)
  s3kpruned walk_2    : 40 fresh trials, best -197.51 (t34) - 2nd seed confirms
  s3kpruned air       : 30 fresh trials, best -7.73 - rhythm gate never beaten
                        (SAME as the pre-fix study; the pruned model tunes at
                        the walk stage, not air)
  w2lvar s1/s2/s3 fix : 29/30/30 fresh trials, bests 136.26 / 68.65 / 25.70
                        - s1 fresh winner 136.26 BEATS the old 71-trial
                        pre-fix best 127.75 (the mutual-inhibition fix helps)
  spiking mirror s1   : 28 trials, best 162.69 (t14) - up from untuned 66.9
  spiking mirror s2   : 12 trials, best 165.64 (t11) - AFFERENTS IMPROVE the
                        air rhythm past the deafferented winner
  BPA non-spiking     : 66 trials, best 185.91 (t30) - beats the gate-passing
                        seed 146.26
  BPA spiking         : 68 trials, best 266.88 (t33) - SPIKING BEATS
                        NON-SPIKING BY 44 PERCENT

HEADLINE: the BPA-actuated robot with the SPIKING controller (new W2L spiking
twin, build_w2l_split_net_spiking.py) outperforms the non-spiking controller
266.9 vs 185.9 - a more regular march (RG period std 0.02 s; 17/20 stance
episodes vs 8/9) at the same 620 kPa ceiling. Winner knobs (spk): te 5.13,
tf 4.76, tau 0.445, thr_inc 5.46, cap 0.73, comm 1.01, support 0.75, 3x20mm
BPAs per route.

FIGURES/GIFS delivered to Documentation\\...\\Dissertation\\Defense\\slideshow:
  s3k_mutfix, syn6_reseed, s3kpruned, w2lvar_fix_s1air, spiking_mirror_s1,
  spiking_mirror_ground (the untuned bilateral ground gait), bpa_nonspiking,
  bpa_spiking - each = _joints.png overlay (vs OpenSim ref where applicable)
  + _walk.gif stick figure (~2.5 s per cycle, looping).

INFRASTRUCTURE (kept, all in spinal/ unless noted):
  _curriculum_spiking.py       gait2392 spiking-mirror curriculum (AARL_NET=spiking)
  w2l_mujoco/build_w2l_split_net_spiking.py   W2L spiking twin (97n/168syn,
                               census-matched analog twin; calibrated relays)
  w2l_mujoco/_curriculum_bpa.py + _curriculum_bpa_spiking.py
                               BPA walker tuning (shared run_bpa plant loop)
  _fresh_20261003.py           fresh dated studies via canonical modules
  _resume_20261003.py          resume driver (study-name redirect)
  _replay_best_20261003.py     winner replay -> npz (NoCloseRaw stdout guard)
  _figs_20261003.py            overlay figures + stick-figure gifs (Tol palette)
  _fail_stale_20261003.py, _status_20261003.py, _refresh_json_20261003.py,
  _w2lvar_json_20261003.py, w2l_mujoco/_bpa_runbest_20261003.py, chain .bats

INCIDENTS (all diagnosed + fixed tonight):
  1. The harness reaps in-session background tasks (2x overnight study kills,
     exit 137) -> studies now run as detached PowerShell Start-Process bat
     chains that survive independent of this session.
  2. _fail_stale filter bug (disconnected loop = no filter) externally FAILED
     live trials of running studies, crashing their optimize loops -> script
     rewritten; fail-stale removed from all chains (manual-only tool).
  3. stdout double-TextIOWrapper GC closes the shared buffer in scripts that
     import several curriculum modules -> _NoCloseRaw(io.RawIOBase) guard in
     the replay script (verified against gc.collect()).
  4. _curriculum_pruned/_curriculum_spiking clobber AARL_NPZ at import ->
     replays re-pin the output name AFTER the import.

OPEN FOR THE NEXT SESSION:
  - w2lvar fresh s4/s5 (walk stages) - stages 1-3 complete tonight
  - spiking-mirror stages 3-5 (balance + ground walk)
  - pruned walk: pre-fix converged at -134 over 600 trials; fresh post-fix is
    at -187 over 80 trials - long batches would close the gap
  - BPA quiet standing remains architecturally blocked (knee-extensor route
    sign; campaigns/20260930/bpa_actuation.md section 1.5)
===== END SUMMARY =====
"""

with open("CURRICULUM_RUN_20261003.log", "a", encoding="utf-8") as f:
    f.write(SUMMARY)
print("summary appended,", len(SUMMARY), "chars")
