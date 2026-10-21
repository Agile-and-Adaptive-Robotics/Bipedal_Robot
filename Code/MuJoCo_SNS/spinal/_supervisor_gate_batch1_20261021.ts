// Supervisor M1 gates for the stage-3 balance + perturbations stream,
// batch 1 of 2 (2026-10-21): stock human model, W2L Hill robot, BPA
// walker - the three models whose results are final. Batch 2 (w2lvar,
// syn6) gates after their s3p studies finish.
//
// The script reads SUPERVISOR_AGENT.md FRESH from the repo and passes
// its full text verbatim as the bulk of each gatekeeper's prompt, with
// a one-paragraph M1 context per the /supervisor-gate procedure.

interface GateItem {
  /** Checklist item from SUPERVISOR_AGENT.md being answered. */
  item: string;
  pass: boolean;
  /** file:line evidence, or the command output that proved it. */
  evidence: string;
}

interface GateVerdict {
  model: string;
  verdict: "SHIP" | "BLOCK";
  items: GateItem[];
  /** Smallest concrete fixes, only when BLOCK. */
  fixes: string[];
  /** Commands the gatekeeper itself re-ran (checklist item 6). */
  commandsRerun: string[];
  notes: string;
}

phase("Read the supervisor contract fresh from the repo");
const contract = await files.read("Code/MuJoCo_SNS/spinal/SUPERVISOR_AGENT.md");

const PY = "D:\\Anaconda\\envs\\myo\\python.exe";
const SPINAL = "D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal";

interface ModelGate {
  key: string;
  label: string;
  claim: string;
}

const models: ModelGate[] = [
  {
    key: "stock",
    label: "human stock s3k model",
    claim:
      "Stage-3 standing balance for the main human (s3k-lineage) model, with perturbations, metabolic cost, and Ib ablation. " +
      "(1) Quiet study curr_s3_balance (optuna_walk.db, log curr_s3bal_20261020.log): 40/40 trials, 40 stand / 0 fall, best = trial 0 (score 17.29) = the incumbent with vestibular gains 0 and rig_scale 1.0 - the searched vestibular gains found NOTHING better (honest negative). " +
      "(2) Fresh-process verify (verify_20261021.log, _s3_verify_stock_20261020.json): sway 0.1095 m, tilt max 28.91 deg, contact sym 0.251, com_z_min 0.802 m, fell=false; recomputed score 17.294510739883755 BIT-EXACT vs DB. " +
      "(3) Rig-share analysis (_s3_rigshare_20261021.json): feet carry 354 N = 50.2% of the 706 N nominal weight, rig holds the other half - a RIG-ASSISTED stand, not free balance. " +
      "(4) Perturbation + metabolic study curr_s3p_perturb (20 trials, 7-eval battery: quiet, 60 N push +-x/+-y, 6 cm platform translate, 6 deg platform tilt): best = trial 0 seed (-46.08), falls 0/7, push-sway mean 0.090 m; the mechanism knobs (stance_width, knee_bias, coact, pelvis_ty) found nothing better than 0 - second honest negative. " +
      "(5) Ib ablation (_s3_ibablate_stock_20261020.json, re-run byte-identical): no fall, but the RIGHT knee drops into 42-50 deg SUSTAINED flexion post-warmup (intact holds +9..+10 deg near extension) and foot load falls 50% to 36% - the Ib-removed model partially crumples into the rig without triggering the fall threshold. " +
      "Claimed verdict: rig-assisted stand with large lean; two honest negatives; Ib removal = unilateral knee collapse-into-rig. Evidence: " + SPINAL + " (curr_s3bal_20261020.log, verify_20261021.log, _s3_verify_stock_20261020.json, curriculum_stage3.json, curriculum_s3p_stock.json, _s3_ibablate_stock_20261020.json, _s3_rigshare_20261021.json, _curriculum_s3p.py, s3p_stock_20261020.log). " +
      "Re-run for item 6 (from " + SPINAL + "): " + PY + " -u _s3_verify_20261020.py stock   (about 70 s; prints the raw bal_* metrics and the recomputed score).",
  },
  {
    key: "w2lhill",
    label: "W2L Hill-muscle robot",
    claim:
      "Stage-3 standing balance for the W2L Hill-muscle robot (w2l_mujoco, bilateralrg2 connectome, M6 gate a). " +
      "(1) Unrigged (rig S=0.0, step regime): FALLS at t=1.18 s (pelvis z<0.6 m path), sway 159 cm, tilt 109 deg - FAIL. " +
      "(2) Latch regime rig 0.0: survives only by hanging - sway_x range 48.1 cm, tilt 4.8 deg, heel duty 0.00/0.00, harness ~+50% weight. " +
      "(3) Rigged (0.25/0.5/1.0, both regimes): all PARTIAL - harness carries ~49-50% of the 411 N weight, heel duty 0.00, tilt 1.0 deg (latch) to 21.4 deg (step). " +
      "(4) Perturbation battery at rig 0.25 (40 N pushes +-x/+-y, 5 cm translate, 5 deg tilt, 12 s): all cases no-fall; post-event sway 2.2 cm (latch) and 2.8-4.8 cm (step). " +
      "(5) Ib ablation (ibnA=0, the Ib group gain into RG-E/InE/PF-E): rigged stands essentially unchanged; unrigged falls at t=1.88 s instead of 1.18 s - falls either way. " +
      "Claimed verdict: the W2L Hill robot does NOT balance; rigged standing is a harness hang with unloaded heels. Evidence: " + SPINAL + "/w2l_mujoco (balance_w2l_20261020.log, balance_w2l_ibablate_20261021.log, _balance_w2l_20261021.cmd, _balance_w2l_ibablate_20261021.cmd, test_w2l_ground.py). " +
      "Re-run for item 6 (from " + SPINAL + "/w2l_mujoco): " + PY + " -u test_w2l_ground.py --phase=stand --dur-stand=12 --sregime=step   (about 30 s; prints the rig ladder incl. the rig-0.0 fall line).",
  },
  {
    key: "bpa",
    label: "BPA-actuated robot",
    claim:
      "Stage-3 standing balance for the BPA-actuated walker (w2l_mujoco/make_bpa_walker.py + test_bpa_stepping.py --balance=1, split-RG W2L drive, 12 BPAs, drop protocol hold 3 s / lower 1 s). " +
      "(1) Support 0.70 (harness carries 70% of 411 N), full 7-case battery (push 40 N +-x/+-y, platform 5 cm translate, 5 deg tilt): every case PASSes the foot-load gate (feet keep >=10% of body weight after the event; measured 10-11%), post-event sway 26-44 cm, root z min 0.934-0.941 m. " +
      "(2) Support 0.85 (push_f only): FAILS the gate at 9%. " +
      "(3) Effort proxy: mean active pressure 272 kPa, peak 310 kPa against the 620 kPa ceiling. " +
      "(4) Ib ablation is N/A: this walker is driven by TONIC open-loop inputs (no afferent wiring exists to ablate) - stated as N/A, not skipped silently. " +
      "(5) The runner itself prints: retained-harness stand, no free-stand claim - this body cannot statically stand, COM sits outside the support polygon. " +
      "Claimed verdict: retained-harness stand that passes the perturbation foot-load gate at 70% support and fails at 85%; no free-stand claim anywhere. Evidence: " + SPINAL + "/w2l_mujoco (balance_w2l_20261020.log, _balance_w2l_20261021.cmd, test_bpa_stepping.py, make_bpa_walker.py). " +
      "Re-run for item 6 (from " + SPINAL + "/w2l_mujoco): " + PY + " -u test_bpa_stepping.py --balance=1 --dur=12 --hold=3 --lower=1 --stand-te=3 --stand-tf=4 --support=0.70 --pushN=40 --pushAxis=+x   (about 25 s; prints the battery block with the PASS line).",
  },
];

phase("Gate each model's stage-3 balance result");
const verdicts = await Promise.all(
  models.map((m) =>
    agent(`gatekeeper-${m.key}`, {
      system:
        "You are a repo-auditing gatekeeper on machine easteregg2 (Windows, Git Bash available). " +
        "You do NOT do the work; you audit it. You must read repo files yourself and re-run at least one " +
        "verification command yourself - never trust the session's claims. Do NOT edit, create, or delete " +
        "ANY file, and do NOT run git commit/push. python is only available as the full path " +
        "D:\\Anaconda\\envs\\myo\\python.exe (bare python is a broken Store stub). If a check is impossible " +
        "to perform, say so plainly in notes instead of guessing.",
    }).ask<GateVerdict>(
      "Milestone: M1. Context: gating the stage-3 standing-balance result for the " +
        m.label +
        ", produced by the 2026-10-20/21 balance campaign in D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal. " +
        "The stream log entries live in reports_20261020/SUPERVISOR_STATUS.md (see the 'stage-3 balance + " +
        "perturbations' entries). The claim under audit:\n\n" +
        m.claim +
        "\n\nAudit it per the contract below. Answer every checklist item with PASS or FAIL plus exact " +
        "file:line or path evidence (items about figures and dissertation .tex are still answerable: check " +
        "that no figure/tex was produced or modified by this campaign - git status and the stream log are " +
        "your friends). Re-run the suggested command (or an equivalent verifier from the repo) and quote " +
        "what it printed. Judge especially whether the claimed verdict row is honest: does the model stand, " +
        "lean on the rig, or game the score?\n\n" +
        "==== SUPERVISOR AGENT CONTRACT (verbatim from the repo) ====\n" +
        contract,
    ),
  ),
);

phase("Collect the verdicts");
const summary = verdicts.map((v) => `${v.model}: ${v.verdict}`).join("; ");
log(`M1 gates batch 1: ${summary}`);
return { batch: 1, verdicts, summary };
