# Easteregg2 history-query prompts (Ben, 2026-10-24, rev 2)

Paste these ONE AT A TIME into the easteregg2 session Ben is using — as of
2026-10-24 that is a FRESH session with NO prior history, so every query
runs against DISK, GIT, and APP INDEXES only, never conversation memory.
If a fact would only exist inside some past conversation, say so and mark
it CONVO-ONLY rather than guessing. (Ben's modeling context is spread
across many conversations on three machines; conversation transcripts are
not assumed recoverable — laptop broken, cloud chats Ben-side, local
ZCode storage holds metadata only.) All outputs go OUTSIDE the repo at
`D:\GitHub\_zcode_backups\machine_context_20261024\` so the working tree
stays clean.

---

## Q1 — Conversion evidence recovery (feeds dissertation SNS §6.3.1)

```
CONTEXT RECOVERY — write to D:\GitHub\_zcode_backups\machine_context_20261024\conversion_evidence.md. Do not commit.

Search this session's history AND this machine's drives/git history for every trace of the Gait2392 / robotbody MuJoCo conversion work: (1) MyoConverter runs on gait2392_robotbody / ConnorBipedal .osim files, especially failures flagged as impossible muscle force/length from BPA routings combined with human muscle parameters; (2) small-muscle pruning (Ben attests 3 per side) — which muscles, which stage, source/output model; (3) patella path handling (moving insertions vs coordinate transform); (4) the OpenSim/ISB frame convention restoration vs MuJoCo's axis convention — the actual transform applied; (5) ground-contact failure symptoms and fixes. Strong leads: Solid_Models\OpenSim\Gait2392_Robotbody\repair_robotbody_muscles.py, Code\MuJoCo_SNS\README.md section "Model repair", Code\MuJoCo_SNS\spinal\DESIGN.md (ISB axes directive + contact diagnostics), and your git log (e.g. commits touching .osim/.py/.xml under those trees).

For EACH finding report: date, exact command or commit, input/output model paths, log or file:line excerpt, and whether it ran on THIS machine (easteregg2) or elsewhere. Tag every row with an evidence level: Ben-attested / code-supported / dynamically reproduced / unresolved / contradicted. Absence is a finding — say so explicitly. Do not infer a fix from the existence of a script.
```

## Q2 — Heavy-run and bridge/report inventory

```
CONTEXT RECOVERY — write to D:\GitHub\_zcode_backups\machine_context_20261024\heavyrun_inventory.md. Do not commit.

Inventory every long computational campaign this session/machine ran or hosts: (1) optuna studies and DBs under Code\MuJoCo_SNS\spinal (optuna_walk.db, optuna_w2lvar.db, optuna_syn6.db, optuna_w2l_contact2.db, others) with dates, study names, winner files, and report folders; (2) Xi/gamultiobj identification runs and Mesh_Optimization sweeps with result mats of record; (3) ANY cross-platform / MuJoCo-Simulink bridge work — CROSS_PLATFORM_REPORT.md, BRIDGE_REPORT.md, E1 reflex / E2 open-loop experiment scripts, W2L matched-circuit records. NOTE: those report files are ABSENT on EB475WS4, and your HEAD~1 commit is "torque paper and bridge" — search your git history and working tree for them explicitly (git log --all --oneline --name-only or similar), and if found, report their full paths and dates.

For each: what ran, when, where the outputs/logs live, and which dissertation numbers (Appendix C Xi values, wean results, torque comparisons) trace to it. Tag evidence level per row. Absence is a finding. Do not rerun anything.
```

## Q3 — Git/branch verification (the "453" mystery is resolved — verify it)

```
CONTEXT RECOVERY — write to D:\GitHub\_zcode_backups\machine_context_20261024\git_forensics.md. Do not commit. READ-ONLY: no pull, merge, rebase, reset, checkout-of-paths, or push (a plain pull is reserved for the campaign session).

RESOLVED CONTEXT (2026-10-24, from EB475WS4): this repo has NO main branch; default is master. This machine is on FEATURE BRANCH KneeTestSetup_BenBo_stw, behind its own upstream origin/KneeTestSetup_BenBo_stw by 4 commits (Ben's 2026-10-23 pushes). The earlier "453 ahead / 6 behind" compared the feature branch against master — normal branch drift, not a sync anomaly. Your job: (1) verify that picture (git remote -v; git branch -vv; git fetch; git rev-list --count --left-right HEAD...origin/KneeTestSetup_BenBo_stw); (2) list the 4 incoming commits (git log --oneline HEAD..origin/KneeTestSetup_BenBo_stw) and what paths they touch — note whether they carry Notes\DISS_STATUS.md, dissertation notes trees, or spinal reports folders this checkout is missing; (3) record the branch layout (which branches exist locally and on origin, their heads and dates) so the campaign session knows the topology; (4) anything else in stashes or other branches that holds dissertation-relevant work. Everything read-only.
```

## Q4R — Cross-conversation residue sweep + session census (replaces Q4)

```
CONTEXT RECOVERY — write to D:\GitHub\_zcode_backups\machine_context_20261024\session_census_and_residue.md. Do not commit.

Ben's neuromechanical modeling context lives in MANY conversations across three machines; transcripts are not recoverable from this fresh session. Recover what those conversations LEFT BEHIND instead:

(1) SESSION CENSUS — enumerate what sessions existed on THIS machine: check C:\Users\<user>\.zcode\v2\tasks-index.sqlite (read-only sqlite; tables tasks / off_peak_tasks carry task_id, title, status, model, created_at, workspace_path) and any C:\Users\<user>\.agents equivalent. Output the full list with dates. Flag every session whose title suggests spinal/W2L/conversion/torque work but which has NO matching artifact found in (2) — those are the conversations whose content may now exist only in Ben's memory or in cloud chats.

(2) RESIDUE SWEEP — walk the durable artifacts those conversations produced in D:\GitHub\Bipedal_Robot: every reports_* folder (all dates), Code\MuJoCo_SNS\spinal\DESIGN.md (its history/gotcha sections), any campaign_INDEX.md / SESSION_NOTES_*.md / *_HISTORY* files, Dig_out logs, tmp\ session scratch. Build a CHRONOLOGY: date → conversation/session (if linkable from census) → ruling or decision made → number produced → file anchor. Tag each row with an evidence level; mark anything that seems consequential but has no disk anchor as UNANCHORED (likely CONVO-ONLY).

(3) End with: the 10 most load-bearing facts that have disk anchors (so the campaign session can rely on them), and a short list of UNANCHORED items Ben may want to re-type from memory or recover from cloud exports.
```

## Ben-only recovery actions (not agent tasks; LOW expectation, Ben 2026-10-24)

- Ben's ruling 2026-10-24: **don't count on the ChatGPT export having much
  value, and don't count on the laptop.** Treat both as zero-yield until
  he says otherwise; the session census + residue sweep (Q4R) is the
  recovery mechanism that matters.
- ChatGPT/Astra conversations are cloud-side; a Settings → Data controls →
  Export data zip CAN be dropped into `D:\GitHub\_zcode_backups\` for an
  agent to mine, but only if Ben ever decides it is worth it.
- Laptop conversations and artifacts wait for drive salvage; nothing
  blocks on them.
