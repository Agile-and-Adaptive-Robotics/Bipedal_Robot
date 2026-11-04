# run_0 termination record (2026-11-04)

## What happened

- Driver `fr_launch_run0` started 2026-11-04 06:41:49 PST (master.log), parpool 10 workers.
- Layer 1 (trunk/deep hip, 5 groups): both branches solved, 10 checkpoints (06:42-06:48).
- Layer 2 (thigh, 10 groups): 5 rigid checkpoints (06:48-06:49), all 10 compliant branch
  solves started 06:49:24-06:49:51 (per-group logs), then **the driver client AND all
  10 pool workers disappeared** (no MATLAB process from the 06:41 launch exists; no
  FAILED lines in any group log). Hard process death ~06:49-06:51, ~2 minutes into the
  layer-2 compliant solves. NOT a graceful stop; the STOP file (created 08:24, after
  death) never reached a live process.
- Preserved work: 15 checkpoints (10 rigid + 5 compliant) = 5 complete groups (psoas,
  iliacus, ercspn, intobl, extobl) + 5 rigid-only (add_mag3, glut_max3, grac, sar,
  vas_int). These are WARM STARTS for the corrected campaign, not results of record.
- Known issues in this formulation (Ben 2026-11-04 directive): constant-Fmax targets
  (no force-length), Xi3 as max-envelope through the Xi0 channel, hook-only CAD
  obstacles, BPA-corridor-as-hard layering constraint, unlabeled aggregate "viol"
  (sar and grac rigid were logged "SOLVED" with viol 1.00e+00 - a 100% violation).
- The layer-2 compliant crash mode is itself a suspect: ALL workers died together
  inside the first compliant evaluations. The corrected campaign MUST smoke-test the
  compliant branch serially in a fresh process before any parfor (fr_xi3_check +
  a one-group compliant solve), and the run-1 driver must not batch 10 compliant
  solves in parallel until that smoke passes.

## Disposition

run_0 is superseded by the corrected campaign (run_1, warm-started from these
checkpoints). STATUS.json was never written (no clean exit). Checkpoint inventory:
`ls run_0/checkpoints` (15 files as of 08:40 PST 2026-11-04).
