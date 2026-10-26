# Dissertation evidence & integration campaign — consolidated easteregg2 prompt (Ben Bolen, 2026-10-24)

This prompt supersedes sending the Revision_D handoff alone: it carries the
machine preamble, the supervision contract, and the independently verified
connectome-audit context, then the handoff directives verbatim.

## 0. Machine identity and hard rules

FIRST run `hostname`. This prompt is written for **easteregg2**. If you are
not easteregg2, stop and read only that machine's AGENTS.md bullet instead.
Facts you will need: repos under `D:\GitHub\`; Python =
`D:\Anaconda\envs\myo\python.exe` (bare `python` is a Store stub; py3.10.21,
mujoco 2.3.7, sns-toolbox 1.5.2, optuna 5.0.0); git CLI on PATH; 10 cores /
128 GB — you are the heavy-compute box; skills at
`C:\Users\Ben\.agents\skills\` (junctions into ZCode_Skills); MiKTeX pdflatex
at `D:\MiKTeX\miktex\bin\x64\pdflatex.exe` (full builds CAN run here).

Rules: the repo is PUBLIC — no credentials, license details, or SSH setup in
anything you write. Ben commits via GitHub Desktop; you do not commit or
push, and you leave the working tree clean (work from copies outside the
repo where directed below). Preserve existing dirty files. Ben is on MOBILE
with a broken laptop: batch your questions, never block on a clarification
you can resolve from the repo, write everything to files, keep chat replies
short with paths. Simscape/smimport and the MuJoCo↔Simulink bridge belong to
the laptop (currently down) — defer that work and note it; do not attempt it
here.

**Reverse SSH is available and authorized** (Ben, 2026-10-24): alias
`eb475ws4` reaches the desktop from here. Use it READ-ONLY: compare live
files, and scp current copies of repo files your checkout lacks into
`D:\GitHub\_zcode_backups\arbitration_inputs\`. Never write to the desktop
except under its `tmp\`.

## 1. Git state — resolved 2026-10-24 (corrects an earlier misread)

The earlier "453 ahead / 6 behind" was a probe artifact: it compared this
machine's FEATURE BRANCH against `master` (there is no `main`; default
branch is `master`). Actual state as of 2026-10-24 ~12:15: this repo is on
branch **`KneeTestSetup_BenBo_stw`** (HEAD beb6c480 "paperwork"), clean
tree, **behind its own upstream `origin/KneeTestSetup_BenBo_stw` by 4**
(Ben pushed those from the desktop on 2026-10-23). No sync anomaly; Ben
fetches/pushes/pulls/commits on both machines.

Do this first: `git pull --ff-only` (safe fast-forward of Ben's own pushed
commits — nothing else). Then re-check the paths §2 flagged missing
(`Notes\DISS_STATUS.md`, dissertation notes trees, spinal reports folders)
— the 4 incoming commits may carry some. Still NO commits, pushes, merges,
rebases, resets, or branch switches by you. The branch-vs-master delta is
normal feature-branch drift; your git history is part of the evidence base
(e.g. HEAD~1 is "torque paper and bridge"). If a referenced file is still
missing after the pull, search history (`git log --all`) and record a
machine-gap, or scp the desktop copy per §0.

## 2. Required first reads (in order; skip-and-record if absent here)

1. The installed **dissertation** skill's SKILL.md and its required
   campaign/status/style reads.
2. Repo-root `AGENTS.md` — dissertation section, Xi semantics, spinal
   contracts, machine guard.
3. `Notes\revision_campaign_20261015\`: BEN_TODO.md, style-contract.md,
   RULINGS_20261016.md, SNS_DISCUSSION_REDESIGN_20261017.md,
   FLOW_REVIEW_20261017.md, CAMPAIGN_REPORT.md.
4. `Notes\DISS_STATUS.md`; newest `Documentation\Project_Digests\digest_*.md`.
5. `CHATGPT_REPORT.md` (the Cxx/Txx/Rxx advisor-thread matrix) — from the
   restored backup in §3, since it is uncommitted on EB475WS4.

Many of these are missing from this checkout by design (see §1); the backup
zip carries the critical ones. Record every gap; continue with what exists.

## 3. Context restore (do this before any work)

Extract `D:\GitHub\_zcode_backups\diss_context_backup_20261024.zip` (44.5 MB,
staged 2026-10-24 by ZCode on EB475WS4) to
`D:\GitHub\_zcode_backups\diss_context_20261024\` — OUTSIDE the repo. Read
its `BACKUP_MANIFEST_20261024.md` first. It contains:

- `Revision_D\` — Ben's prose working copy of record (changed chapters:
  Introduction, Background, Database, SNS; `review\` holds before-copies,
  `changes.json`, `changes.diff`; `ZCODE_HANDOFF.md` is the handoff text).
  Uncommitted on EB475WS4; treat as the newest prose source. `thesis.bib`
  traveled binary-safe in the zip — never round-trip it through
  default-encoding PowerShell; Latin-1 byte-safe edits only. `build_*` dirs
  were excluded on purpose: rebuild fresh.
- `CHATGPT_REPORT.md` — the live matrix you will update.
- `connectome_audit\` — AUDIT.md, ZCODE_HANDOFF.md, inventory.json.

If `D:\GitHub\_zcode_backups\machine_context_20261024\` exists (Ben ran the
2026-10-24 history queries), read those files first — they carry this
machine's conversion evidence, run inventory, git forensics, and the
session census + residue chronology (Q4R). Context reality (Ben,
2026-10-24): his modeling work is spread across MANY conversations on all
three machines; conversation transcripts are NOT assumed recoverable
(fresh query sessions; laptop broken; cloud chats are Ben-side). The
durable record is disk artifacts + git — treat anything without a file or
commit anchor as UNANCHORED and surface it to Ben rather than using it.

Your own working outputs for this campaign go to
`D:\GitHub\_zcode_backups\diss_campaign_20261024\` (create it). Nothing in
the repo tree gets edited except where PART B explicitly directs prose edits
in the restored Revision_D copy.

## 4. Supervision contract (supervisors are IN the workflow)

Work in batches: **arbitration** (§5), then **BW** the W2L connectome fix
(§5's Batch BW — Ben, 2026-10-24); then **B1** conversion evidence (PART B
Priority 1 + 1.1);
**B2** five-system status + existing-results audit (Priority 1 item 2–3);
**B3** one exact connectome + published protocol (Priority 2); **B4**
database checks; **B5** editorial small items. After each batch's rows are
drafted: run the read-only **diss-supervisor** agent if installed on this
machine — under **THREE gates only** (Ben, 2026-10-24: external plagiarism,
AI-slop patterns, factual drift, comment-intent — **the self-plagiarism
gate is REMOVED**; a stapled-articles dissertation is acceptable here and
the torque paper may already be a pre-print by acceptance, so text
overlapping Ben's own published/submitted papers is NOT a finding). If the
installed agent's gate list is fixed and includes self-plagiarism, run it
but mark any self-plagiarism findings as WAIVED per this ruling, recording
the ruling; otherwise run the three gates manually per the dissertation
skill and record the check. A matrix row closes only on supervisor GO plus
the evidence required by PART B's Return contract. Repeated concern reopens
the row. No GO from a clean build or a populated checklist alone. Preserve
old scores and their evidence in CHATGPT_REPORT.md; no silent deletions or
number changes anywhere.

Standing dissertation rules that hold here too: every number traceable to a
Testing_Data file, generating script/mat, or citation; per-pass verification
build uses the full pdflatex path, a fresh `build_*` dir, and bibtex with
absolute BIBINPUTS/BSTINPUTS; no machine-voiced comments in .tex; no dates
or campaign provenance in prose; Modern Robotics notation; engineering
diameter symbol (mathabx `\diameter`), never `\phi`; "faster than real time"
only in its scoped, literature-supported form; minimal AI self-narration in
repo docs; final prose polish is ChatGPT Astra's (bounded, evidence-checked
sections) — GLM-5.3/5.3-Flash (this ZCode setup) does testing, modeling,
optimization.

## 5. Verified connectome-audit context (ZCode on EB475WS4, 2026-10-24)

The ChatGPT connectome audit in the backup (`connectome_audit\`) was
independently machine-verified on EB475WS4. Treat it as accurate WITH these
amendments:

- AUDIT.md finding 2's phrase "four phase cells and six NMF synergies"
  MISATTRIBUTES syn6's six NMF channels to the Di Russo draft. The draft
  (`replication\dirusso_connectome_draft.json`) actually specifies a
  two-phase oscillator pair (phi_L/phi_R) and FIVE raised-cosine primitives.
  Apply this correction anywhere the audit is quoted; never propagate the
  error into the dissertation.
- Verified source conflicts of record: Di Russo RC polarity — Methods say
  antagonist-RC mutual inhibition, Table 1 caption says reciprocal
  excitation; the repo keeps mutual inhibition, flagged as a real source
  inconsistency. Shinohara Table A.1 (PDF p.31) prints RG-E→PF-F at 0.7
  with a dash under PF-E; Ben ruled RG-E→PF-E with the printed-table
  conflict retained in provenance, tested only as a labeled sensitivity.
- Shinohara's tau_h(V) = 320 + 320/cosh((V+35)/15) ms must not be replaced
  by the fixed-tau workaround; a bell-tau SNS backend already exists in
  `replication\shevtsova_sns.py`.
- The corrected W2L template (`bilateralrg2`, census 96/204/12) is missing
  from the editor HTML's embedded library and the generator's menu — fixed
  by Batch BW below (Ben, 2026-10-24: stop deferring it).
- The syn6 synergy-timing challenge (Ben's open objection, preserved in the
  audit's ZCODE_HANDOFF.md with exact verified numbers: W 43×6 per side;
  families S1 F, S2 E, S3 E, S4 F, S5 F, S6 E on both sides; bilateral-mean
  stance fractions 0.195713/0.952928/0.922159/0.246330/0.410065/0.889709;
  S5 mixed-band exception with default-off reciprocal inhibition) is
  unresolved. It is NOT validated physiology; never present it as such.
- PART B Priority 2 candidate selection: on-disk execution evidence exists
  for Shevtsova (`replication\verify_shevtsova.py`, `verify_shevtsova_sns.py`)
  and Rybak (`verify_rybak2006.py` + `verify_rybak2006.log`). Check the
  logs' dates and content FIRST — the older README caveat says these
  translations were not built/verified as running models; resolve which is
  current before choosing. Note the reports_20261016/20261020 folders may
  be absent here; read them if present, else record machine-gaps.
- Spinal contracts in any verification: default census + the bit-exact
  anchor −161.56754173676563 (s3k stage-4, full_rules=0); fresh-process
  verification for any wean/basin claim; MuJoCo 2.3.7/implicitfast/2 ms
  frozen through the dissertation; frozen/NaN runs never count as walking;
  no harness-free walking claim anywhere in the program's records.
- Optimization record cautions that are already repo-documented: the
  2026-10-20 corrected-W2L drive ladder is withdrawn (dscale was
  bookkeeping; te/tf stayed free), best strict basin 9/12 vs required
  11/12 (marginal phase robustness, NOT robust walking), and the
  AnimatLab tuned-winner transfer fell and latched (no drop-then-walk).
- CROSS_PLATFORM_REPORT.md, BRIDGE_REPORT.md, and the E1/E2 scripts are
  ABSENT on EB475WS4 — search THIS machine's tree and its git history
  (HEAD~1 is "torque paper and bridge") before concluding anything
  about them.

### Arbitration task: ZCode's verification vs ChatGPT's audit (Ben, 2026-10-24)

Ben asked THIS machine to independently check who is right. Re-verify each
disputed point with your own commands (files: `connectome_block_editor.html`,
`connectome_templates.json`, `make_editor_templates.py`,
`replication\dirusso_connectome_draft.json`, `replication\shinohara_rules.json`,
`build_network_syn6.py`; the audit itself is in the backup zip). FIRST check
freshness: `git log -1 --format=%%ci -- <path>` on each; if a file is
missing or older than 2026-10-06, scp the desktop's current copy via
reverse SSH into `D:\GitHub\_zcode_backups\arbitration_inputs\` and arbitrate
against that. Return a verdict per point — ZCode correct / ChatGPT correct
/ both partially correct / indeterminate — with the command and output:

1. AUDIT.md finding 2's "four phase cells and six NMF synergies" — ZCode
   claims the Di Russo draft actually specifies a two-phase oscillator pair
   and FIVE primitives, and that the six NMF channels belong to syn6.
2. `bilateralrg2` absent from the editor HTML's embedded library and the
   generator's menu while present in the sidecar (24 embedded entries, all
   matching the sidecar).
3. Shinohara Table A.1 (PDF p.31) prints RG-E→PF-F at 0.7, dash under PF-E;
   Ben's ruling is RG-E→PF-E with the conflict retained.
4. Di Russo RC polarity: Methods say antagonist-RC mutual inhibition; Table
   1 caption says reciprocal excitation — a real internal source conflict.
5. syn6 numbers: W 43×6 per side; families S1 F, S2 E, S3 E, S4 F, S5 F,
   S6 E both sides; bilateral-mean stance fractions
   0.195713/0.952928/0.922159/0.246330/0.410065/0.889709.
6. The withdrawn wean ladder, 9/12 best strict basin vs required 11/12, and
   the AnimatLab tuned-transfer fall/latch — confirm these match the
   reports of record (if those folders are absent even after the §1 pull,
   say so rather than confirming).

Write the verdicts to `D:\GitHub\_zcode_backups\diss_campaign_20261024\audit_arbitration.md`
and lead your first status report with them so Ben can close the question.

### Batch BW — fix the W2L connectome (Ben, 2026-10-24: "I DO want it fixed. I'm tired of that being put off.")

Run right after the arbitration, before B1. All edits happen in the repo
working tree ONLY where Ben's originals are untouched; work on copies for
anything inside the .aproj. AnimatLab headless is PROVEN on this machine
(the 2026-10-20 drop report ran here).

**Tier 1 — editor sync (mechanical, just do it):**
1. Add `bilateralrg2` to `make_editor_templates.py`'s menu mapping and
   regenerate/synchronize so the editor's EMBEDDED library carries it
   alongside the sidecar (all shared entries must remain byte-equal;
   census 96/204/12 must match `spinal\w2l_cpg\`'s miner).
2. Verify: the editor presents `bilateralrg2` (embedded library wins over
   sidecar on load — that is the whole bug); no duplicate labels or
   dangling endpoints; regression-check one other template renders
   unchanged. Per AGENTS.md, the editor's localStorage is per-origin:
   clean up any TEST_* artifacts you create, and do not clobber Ben's
   saved canvases.
3. `smoke_w2l.py --template bilateralrg2` must still PASS, and the old
   `bilateralrg` regression must remain bit-identical (1.000 s / −0.525 /
   22-22 per reports_20261020\w2l_contact2\README.md).

**Tier 2 — solve the MuJoCo wiring: Optuna over conductance classes
(Ben, 2026-10-24: "shouldn't the Optuna optimizer solve for all
conductances?" — yes):**
4. Target = the SNS/MuJoCo W2L network (the `bilateralrg2` template →
   `spinal\w2l_mujoco\test_w2l_ground.py --net=bilateralrg2`), NOT the
   .aproj. First, list the template's documented source-sign and
   conductance-mapping deviations and close each: either it carries a Ben
   ruling (F1 transcribes INHIBITION per his 2026-10-03 export — the sign
   stays) or it is a transcription bug, fixed in the template now.
5. Parameterize conductances PER TARGET — per-synapse where the pathway
   demands it: every MN input conductance and every afferent feedback
   strength gets its own value (Ben, 2026-10-24: family-class collapsing
   of MN inputs and afferent feedback into shared strengths is
   OVERRIDDEN). Seed literature values where they exist (Shinohara
   Table B.3 per-MN PF weights; Di Russo signed per-muscle W; eq 10/11
   group gains) and let Optuna solve the rest as individual bounded
   dims. The anti-lottery machinery (basin gate, gain-sensitivity probe,
   remove-and-retest) applies as VALIDATION of the per-conductance
   winner — it is NOT a reason to shrink the search space. Signs and
   topology remain fixed per Ben's rulings. Follow the JSON RULE: every
   new params.G knob is saved in the winner json with an `if key in
   best` loader branch, and default values keep the current network
   bit-identical. Mining/templates must preserve per-instance values with
   the correct KIND: per-synapse SynAmp (µS) for chemical synapses;
   per-path amplitude (nA) for afferent current injections; MN→muscle
   couplings are activation mappings, not synapses. The 2026-09
   "effective-strength doctrine" (type-level conductance only,
   per-instance = GUI prefill) is OVERRIDDEN by this ruling.
6. Run the study on this machine (db in the spinal optuna family, e.g.
   `optuna_w2l_condsolve`; metric + survival gate VERBATIM from
   reports_20261016\w2l_wean / reports_20261020\w2l_contact2). PIN the
   drive (Ben's literal te/tf defaults or one stated fixed value) so the
   optimizer cannot compensate through drive — that free-drive flaw is
   what invalidated the Oct-20 ladder. No wean-ladder claims in this
   batch.
7. Acceptance = the strict basin gate (12 fresh-process trials, ±1% on all
   knobs, pass bar 11/12) plus inter-leg phase metrics — the documented
   fragile piece (best so far 9/12, with the ch_contra 0.5–0.8 plateau).
   On finalists: per-term diagnostics, gain-sensitivity probe,
   remove-and-retest for lottery weights, and FRESH-process verification.
   Frozen/NaN runs fail. No harness-free walking claim regardless of
   outcome. If the solve plateaus below 11/12, report it as DESIGN.md
   open-issue-#4 evidence for a restructure decision — do not silently
   accept a marginal winner.
8. The .aproj's own F1/V3-type discrepancies stay flagged for Ben's GUI
   pass someday; they are out of scope here.

Deliverable: editor fixed + verified (Tier 1); deviation table closed;
the conductance-solve study with winner json, strict-basin outcome,
anti-lottery checks, and an honest verdict; the s3k default census and
bit-exact anchor −161.56754173676563 re-verified untouched; Ben's
originals byte-untouched. Surface any choice needing Ben as one-word
mobile options after applying documented defaults.

## 6. PART B — handoff directives, VERBATIM

# ZCode handoff: evidence and integration work after the local prose pass

Read CHATGPT_REPORT.md first. Ben authorized local editing and asked for this handoff; Overleaf remains read-only. Do not upload or replace chapters from this folder without a new live-source comparison. Do not commit or push. This working copy is Revision_D; its changed chapters are Introduction, Background, Database, and SNS. Use review/changes.json and changes.diff to merge targeted edits, not whole files. Frozen before copies are in review/before.

## Start here on the other machine

This handoff is self-contained. Paths below are repository-relative; do not assume the desktop D: path exists on the machine running ZCode. FIRST run hostname and read only that machine's AGENTS.md bullet. Locate that machine's Bipedal_Robot and ZCode_Skills repositories. Use its documented Python/MATLAB environments. Preserve existing dirty files. The prose revisions and report are local/uncommitted on EB475WS4: if Revision_D is absent, obtain Ben's synchronized files before claiming to have edited or built them. Do not overwrite another machine's newer source.

Load and follow the installed **dissertation** skill first, including its required campaign/status/style reads. Use **mujoco** for conversion, actuator, joint/frame, and contact verification; **animatlab** for the original Li/Deng/W2L models and reference experiments; **matlab** for BPA equations and Simulink/Simscape work; and **neuroscience** for source-paper circuits, physiological interpretation, and replication protocols. Use the **sns-toolbox** skill if available for network construction and implementation comparisons. Discover actual skill locations and read their SKILL.md files; do not merely say the names in a plan. Use relevant installed plugins/connectors for source documents, literature, repository evidence, and model tools, subject to access and Ben's existing authorization. If a requested skill/plugin is unavailable, record which task is affected and use a supported alternative; do not invent tool access or claim execution on an inaccessible machine.

Ben's approved model routing: GLM-5.3 / GLM-5.3-Flash for testing, modeling, and optimization in his existing ZCode setup; ChatGPT Astra for the later bounded prose pass. This document does not authorize inventing results or optimizing around transcription errors. No credential setup is needed just to read this handoff.

## Routing and acceptance

Ben wants GLM-5.3 / GLM-5.3-Flash for testing, modeling, and optimization through his ZCode setup. Use the existing configured access, not a newly created provider integration. Final prose polish is for ChatGPT Astra, limited to evidence-checked changed sections. Do not spend a premium-model pass on broad repository rediscovery or loop indefinitely polishing the same paragraph. Do not optimize before verifying transcription and baseline behavior.

## Priority 1: provide actual findings, not more abstract prose

1. **Conversion checks (SNS §6.3.1).** Inspect repair_robotbody_muscles.py and the generated/source model comparisons. List the actual muscle-path/joint issues corrected, if any, and the exact checks that passed. Return source file/function, before/after model property, a rerunnable check, output, and proposed two-sentence result. The advisor's wrong-origin/left-hand-axis examples were hypothetical, not findings. Do not label a procedure successful from its existence alone.
2. **Five-system status (SNS §6.2 / Table 6.1).** For each of the five tasks, give implemented model, neural controller, mechanical model, actuation, input protocol, output evidence, and explicit remaining gap. Confirm bilateral PF organization and contributors for the old 'for each leg?'/'We?' comments. The local pass retained the captured source's I/Morrow attribution; it did not independently establish authorship history. The fifth task is not complete merely because the SolidWorks mechanism imports. Check the rigid ankle limitation and whether neural/BPA/mechanical coupling actually ran.
3. **Existing results first.** Audit CROSS_PLATFORM_REPORT.md and its referenced logs, BRIDGE_REPORT.md, E1 reflex and E2 open-loop experiment scripts, and W2L matched-circuit records. Separate open-loop 92-actuator drive, closed-loop reflex exchange, driven rhythm-core agreement, and full neuromechanical replication. Report which actual source/protocol matches each claimed result. Do not rerun tests already adequately evidenced unless a specific mismatch requires it.

## Priority 1.1 supplement: Ben's actual conversion history and replies

**Do not treat the problems as merely hypothetical.** The advisor's illustrative wrong-origin/left-hand examples were hypothetical; Ben's following replies describe actual problems. The previous prose pass had not captured these replies and consequently omitted important specific context. Read them together before revising SNS §6.3.1.

**Advisor, 23 October, 11:31 a.m.** Asked what problems the conversion checks actually found and corrected, rather than implying defects by calling the procedure an "audit."

**Ben's Overleaf reply, 23 October, 2:03 p.m. (verbatim):**

> I ran into problems but I wasn't sure how much to detail here. Some of the very small muscles had to be pruned (3 per side), Mujoco uses the "wrong" frame axis convention compared to OpenSim and the International Society of Biomechanics (ISB). The patella being defined by muscle insertion locations as a function of knee angle in the femur frame instead of via cordinate transform caused problems.

**Ben's Overleaf reply, 23 October, 2:05 p.m. (verbatim):**

> Oh yes, and I first shot for the moon by importing the gait2392_robotbody from my master's project (BPA routing but human muscle parameters) which led to errors due to Myo's checks against impossible muscle force.

**Advisor, 23 October, 11:32 a.m.:** "What BPA interface? Non seq"

**Ben's Overleaf reply, 23 October, 2:07 p.m. (verbatim):**

> I need to clarify that I tried to import the robot (BPA routings on human frame) first or at least in addition to bringing in the gait2392 human biomechanical model.

**Ben's additional instruction in this chat:** the model was changed back to the OpenSim/ISB coordinate system rather than the MuJoCo convention; ground contact caused an issue; initial robotbody import combined BPA routings with unmodified human muscle characteristics and triggered impossible-configuration flags; muscle repair and other corrections were performed. Ben is uncertain whether particular work ran on EB475WS4 or easteregg2. His chat said "ISM"; his explicit Overleaf reply identifies **ISB**, the International Society of Biomechanics. Use the explicit expansion, and verify the actual frame transform.

### Evidence recovery tasks and required result

| Issue confirmed by Ben | Recover from code/logs | What to write after checking |
|---|---|---|
| OpenSim/ISB frame restored | Actual transform, axis ordering, up/forward axes, handedness, frame of meshes/joints/forces, gravity/contact-plane consistency, and before/after comparison. Distinguish a global frame change from a joint-sign change. | State exactly what convention was retained/restored and how preservation of intended motion/torque was checked. Do not equate a different axis convention with left-handed coordinates. |
| Ground-contact issue | Failure symptom, model/version, collision shapes/meshes, foot geometry, gravity/ground orientation, solver/contact settings, corrective change, and rerun. Search this separately from later controller tuning. | State actual cause if demonstrated; otherwise state symptom and remaining uncertainty. Do not assume the frame change caused or fixed contact without evidence. |
| Robotbody import before/in addition to human reference | Exact input .osim, date/run/machine, which MyoConverter stage/check failed, error text, path/length/force assumptions, and whether model or muscle parameters were altered. | Explain that BPA routing with human muscle characteristics was the initial incompatible combination; separate it from conversion of the human Gait2392 reference model and later custom BPA force coupling. |
| Small-muscle pruning: 3 per side in Ben's reply | Names of removed muscles, precise source/output model and stage, rationale, and actuator-count reconciliation. | Keep the source, intermediate, and final counts separate. Do not silently revise the production 92-actuator claim or assume these are the six removed actuators in every model. Flag disagreement for Ben. |
| Patella path definition | Moving insertion locations as functions of knee angle in the femur frame; the actual representation expected by converter; any explicit patella body/joint or path repair introduced; geometry/kinematics checks. | Name what failed and the implemented remedy. Do not claim that a moving insertion is itself biologically wrong or assert that a particular repair solved it without a trace. |
| Other muscle repair | Actual repair script location, source/output model hashes or identifiers, origin/insertion/via-point changes, and geometry/torque-sign checks. | Give one or two concrete verified repairs and their outcome, not a list of vague "audit steps." |

First leads: Code/MuJoCo_SNS/README.md §2 "Model repair — gait2392_robot.osim"; search the repository for repair_robotbody_muscles.py instead of assuming it lives directly under Code/MuJoCo_SNS. The repair script is at Solid_Models/OpenSim/Gait2392_Robotbody/repair_robotbody_muscles.py. The README reports mirroring the intentionally modified right-side paths and re-enabling test-disabled muscles. Verify the actual script and output, and distinguish these repairs from the small-muscle pruning and impossible-force checks in the replies. Code/MuJoCo_SNS/spinal/DESIGN.md contains Ben's ISB/OpenSim axes directive and related contact diagnostics. Check git history, reports, machine-local scratch logs, and the generated model tree. Use relevant plugins and the named skills. Respect AGENTS.md machine guards; do not assume absence on one machine means the work never happened.

Return a table with **source model → failure/check → fix → verification → residual limit → evidence path and machine**. Mark every item as Ben-attested, code-supported, dynamically reproduced, unresolved, or contradicted as appropriate. These are distinct evidence levels. Existing successful logs may suffice; rerun only where necessary. Then propose a short conversion-results paragraph that explains what the checks found, preserves the distinction between initial robotbody and human-reference conversion, and ties the later BPA-force interface to that sequence. Keep extended implementation detail in the appendix/report. Do not omit Ben's context merely because the final fix is still being located.

## Priority 2: one exact connectome and one published protocol

Read Code/MuJoCo_SNS/spinal/DESIGN.md before touching runner.py/build_network.py/params.py. Preserve the default census and bit-exact anchor contract; use fresh processes for relevant verification. Keep MuJoCo 2.3.7 and the recorded step through this work. Heavy runs belong on easteregg2; Simscape runtime work belongs to the laptop.

Start by auditing replication/README.md, the selected source paper and supplement, and the actual current implementation. Do not take the README's dated draft status or the presence of verify_*.py as the final execution verdict. Check newer logs/artifacts on the machines ZCode can access.

For ONE selected source model, produce:

- Node table and directed edge table with source, target, sign, weight, convention, delay, cell type, and exact page/table/figure evidence. Keep a separate unresolved-edge list. Do not infer an obscured connection from what seems physiologically reasonable.
- Figure crops or overlays for disputed edges, cross-checked against supplementary tables/equations. Compare the built SNS network programmatically with the approved edge table: missing edges, extra edges, sign errors, and weight differences, not just counts.
- Baseline traces under the published drive, initial conditions, mechanics, and solver protocol. Then ONE specified stimulus or deletion experiment, with paper target, measurable output, justified tolerances, and raw traces. An oscillating isolated RG is not evidence for a coupled walking replication. A port on different mechanics is an adaptation and must be labeled that way.
- Only after topology and baseline pass: justify a bounded parameter search if needed. Do not tune around an uncertain sign or omitted pathway. Frozen/NaN runs fail. No harness-free walking is established in the program's records.

Use the existing order in replication/README.md as context, but select the best evidenced feasible next source rather than promising all replications in one pass. Return a concrete verdict and a blocker if the source specification cannot be recovered.

## Database checks

The local prose pass verified citation direction in build_citation_graph.py and build_app.py and app-only archive status in sadb_app.html. Verify actual live-base schema and current export coverage before asserting complete pathway/paper/review relationships. Supply one real linked example, with record IDs and source evidence, if permitted. Do not fabricate an example relationship or claim a gap is filled just because the record types exist. Keep recorder records and all corpus values intact.

## Small unresolved editorial items

- Exact highlights for 13 collapsed older anchors are not available; recover them from review history or use the explicitly labeled historical highlights. Do not treat EOF as the original section.
- Background's remaining ethical-case and scientific-validation claims need a scoped fact check: broad replacement of animal experiments, unsupported validation-chain claims, and claims that all interventions become routine should not be polished into stronger assertions. The present pass fixed repeatability and throughput wording only. There is a non seq. that leads from one paragraph to the cruelty described by the figure in the next paragraph. Bridge the gap. Ben has been vegan for over half his life. He's not afraid to fight his advisor on the subject of animal rights.
- Background Haeufle robustness wording: check the actual paper before replacing a technical result with a plausible paraphrase.
- Earlier-library comparison: identify the actual Hunt/Fletcher libraries meant by the advisor before comparing them. Do not substitute a different library or remove a valid requested comparison because the citation is missing. Edit: We don't have these simulink libraries. Idk what he's talking about. This is a good example of you hallucinating an answer. Learn from this.
- Build support outside the five captured chapters came from Revision_C. Reconcile the unedited chapters, main.tex, bibliography, and figures with a fresh complete Overleaf export before any upload. Bibliography was copied byte-for-byte, not re-encoded. The existing Flx20mm_meetsHuman.pdf had to be copied from tmp/overleaf_live_20261023 into this working copy for compilation.
- The build has 26 overfull boxes, predominantly outside this targeted prose work; locate and disposition them. Do not claim clean layout from zero compilation errors.

## Return contract

For every handled matrix row, return: original request, section/line and exact before/after, evidence and command/log, result, remaining limit, and self-score. A clarity revision can earn 3/4 after these checks; 4/4 requires Ben/advisor acceptance. Repeated concern reopens the row. No GO based solely on a build or a populated checklist. Update the matrix and CHATGPT_REPORT.md, preserving old scores and their evidence. No administrative report text belongs in the dissertation.

## 6a. Ben's amendments, 2026-10-24 (SUPERSEDE conflicting PART B lines)

PART B above stays verbatim as the record. Where these amendments conflict,
the amendments win:

1. **Self-plagiarism gate removed** from supervision (see §4). Text
   overlapping Ben's own published/submitted papers is not a finding;
   a stapled-articles dissertation is acceptable, and the torque paper may
   be a pre-print by acceptance.
2. **Ethical subsection of Background: NO revision.** Skip PART B's
   "Small unresolved editorial items" bullet about the ethical-case /
   animal-rights bridging entirely — leave that text untouched. The OTHER
   editorial bullets (collapsed anchors, Haeufle wording, library
   comparison = ask-the-advisor, Revision_C reconciliation, overfull
   dispositions) still stand.
3. **Table 6.1 is not the metric.** Ben is confused by Table 6.1 itself, so
   B2 must BEGIN with a plain-language explanation of what the table
   currently says, what its five rows/tasks are, and exactly what is
   ambiguous — BEFORE any completion work — and must disambiguate Ben's
   question "is anything here the five-part curriculum?" Candidates to
   check and distinguish: (a) the five tasks of SNS §6.2 / Table 6.1;
   (b) the five-STAGE optuna curriculum family in the spinal stack
   (curriculum_w2lvar_stage5.json, curriculum_syn6_stage5.json winners);
   (c) anything in SNS_DISCUSSION_REDESIGN_20261017.md. Present the mapping
   in one short table for Ben's mobile review; make NO prose edits to
   §6.2/Table 6.1 until he confirms what the five things are.
4. **Arbitration first:** run §5's ZCode-vs-ChatGPT arbitration before B1.
5. **Reverse SSH to eb475WS4 is authorized** (§0) — use it for freshness
   checks and read-only copies instead of declaring machine-gaps early.
6. **The W2L connectome WILL be fixed in this campaign** (Batch BW, §5) —
   on the **MuJoCo/SNS side**: the editor gets `bilateralrg2`, and the
   corrected W2L net's conductances are SOLVED by an Optuna conductance-
   class study on this machine (signs and topology stay per Ben's rulings;
   only magnitudes are searched). Ben's .aproj originals are untouched —
   their F1/V3-type discrepancies remain flagged for his own GUI pass.
   This supersedes §5's earlier line deferring the `bilateralrg2` editor
   gap.
7. **Simscape/Simulink capability reality (dynamically verified
   2026-10-24 on EB475WS4; report:
   `Code\Matlab\SNS_Simscape\logs\desktop_capability_20261024_REPORT.md`).**
   The base-Simulink SNS stack is VERIFIED on EB475WS4: library regen
   R2025a-native, actuator tests 4.9e-10 N, units test 6.07e-04 mV, all
   three demos at their documented numbers in fresh processes, the
   410-neuron production network verify 4.2e-06 mV, winner slice all 4
   checks. Simscape Multibody checkout is BLOCKED on this machine AND on
   EB475WS4 (re-verified today by `sns_urdf_smoke.m`: smimport dies at
   the license callback) — the smimport routes stay laptop-only until Ben
   gets a better MATLAB license somewhere; do NOT attempt them. Known
   open bugs: the Deng CPG BUILDER fails on R2025a (`pair1_s1/2`
   generated-name error; it deletes-then-fails) AND the committed
   SNS_Deng_CPGDemo trio does not compile on R2025a (2-wide/1-wide RG→PF
   port mismatch inside the demo; RG/PF themselves are correct 1-wide —
   probe: `Code\Matlab\SNS_Simscape\dev\probe_deng_ports_20261024.m`;
   the September 1.938 s verification stands as a laptop-era record;
   repair is a bounded campaign task). The winner-slice builder mis-paths
   `SNS_Library_R2025a.slx` (its RUN verification passes). MinGW64 is
   configured on EB475WS4, so rebuilding the MuJoCo↔Simulink bridge mex
   there is feasible and is the sanctioned next Simscape-side step.
   Simscape work you CAN do here: anything base-Simulink — build models,
   verify nets, export figures — one model per FRESH process (in-session
   runs right after a library rebuild fail).
8. **THE FIVE-PART CURRICULUM (Ben's definition, 2026-10-24 — this IS the
   "five-part curriculum" asked about; it is NOT Table 6.1's five tasks
   and NOT the optuna stage-5 family; B2 must present the distinction).**
   Stage 1 and 2: AIR STEPPING. Stage 3: BALANCE, no walking — 3a)
   balance with a harness, 3b) balance without harness after weaning.
   Stage 4: WALKING WITH A HARNESS. Stage 5: walking without a harness.
   The harness support must MIMIC THE LITERATURE — depending on the
   model, that often means constraining to sagittal-plane locomotion and
   body-weight-support percentages reduced across stages (clinical
   Lokomat/BWSTT practice: ~30–50% initial, progressively reduced — see
   `Code\Matlab\SNS_Simscape\reports_20261024\literature_values_reflex_curriculum.md`).
   Ben's order on the air-stepping/locomotion work: "Both. All." — MuJoCo
   AND Simulink. Mapping of existing records: W2L +1.718 m/30 s
   harness-supported antiphase stepping = stage 4; the human-model
   in-place stepping weaned to drive 0.0 is NOT stage 5; no stage-5
   result exists — say so. Stages 1–3 have not been run under this
   curriculum and are new work (Batch BC).
9. **Demo validation work order (Ben, 2026-10-24: "Run the SNS toolbox
   tutorials and the AnimatLab tutorial and match their responses. Do a
   literature search and find values from actual studies.").**
   Tutorials 1/2/4/8 were re-run fresh on EB475WS4 2026-10-24
   (`spinal\sns_tutorials\rerun_no_render_20261024.log`; render skipped —
   no graphviz there). Still to do: numeric response extraction + Simulink
   SNS twins of the tutorial circuits; obtain and headless-run the
   AnimatLab Human Stretch Reflex tutorial
   (animatlab.com/Help/Tutorials/Examples/Human-Stretch-Reflex; ModelDB
   Szczecinski simple-stretch-reflex) on EB475WS4 and match responses;
   verify the literature values at their full texts before any citation.
   First real comparison already in hand: beer demo reflex ON/OFF gives
   ~3.9x apparent stiffness (reflex share ~74% of total) vs literature
   30–60% for human ankle — same order, same mechanism; tune the demo's
   series compliance to human values before stronger claims.
10. **Deng + Li rebuild (Ben, 2026-10-24: "THEN FIX IT!! You can
   literally rebuild the Deng and Li models from AnimatLab into SNS
   Simscape. They're crude blocks ffs.").** Rebuild BOTH as clean
   base-Simulink models from the AnimatLab sources: Li from the
   `li` template mine (58n/98e, `spinal\connectome_templates.json`);
   Deng from the Nourse tables + Biped_2xCPG audit. Replace the broken
   R2025a Deng trio. Verify each against the AnimatLab originals run
   headless (AnimatSimulator on EB475WS4) and the numpy Deng reference
   (period 1.938 s).
11. **Dissertation demo images (Ben, 2026-10-24: the current images are
   "terrible").** Redraw the KneeReflex/CPG leg animation to look like an
   actual leg and the beer demo like an actual arm holding an actual cup;
   replace the images used in the dissertation (60-sns components) and
   regenerate under the campaign style contract.
12. **NO FAMILY-CLASS CONDUCTANCE COLLAPSING (Ben, 2026-10-24 — issued
   after catching this in progress; supersedes anything above that
   implies it).** MN input conductances and afferent feedback strengths
   are PER-TARGET: per-muscle, per-path, individual values. Every model
   built or tuned in this campaign carries individual conductances —
   seeded from literature where values exist (Shinohara B.3, Di Russo W,
   eq 10/11), Optuna-solved otherwise, never copied from a sibling
   because they share a family label. Inventory and FLAG (do not delete)
   every existing build that uses shared family strengths, for per-target
   rebuild. Template mines must record per-instance values with the
   correct KIND (per-synapse SynAmp in µS; per-path nA for afferent
   current injections; MN→muscle couplings are activation mappings, not
   synapses); the 2026-09 effective-strength doctrine is overridden. The
   Deng/Li
   rebuilds (amendment 10) carry the AnimatLab originals' per-synapse
   values (Li: 1–8 µS per synapse).

## 7. Deliverables checklist (dissertation-facing; update as you go)

1. Arbitration verdicts: ZCode-vs-ChatGPT, per-point with commands.
2. BW: W2L connectome fixed and verified on the MuJoCo side — editor
   presents bilateralrg2 (smoke + bit-identical regressions); deviation
   table closed; Optuna conductance-class study with winner json,
   strict-basin outcome (11/12 bar), anti-lottery checks, fresh-process
   verification; honest verdict; s3k anchor untouched; originals
   untouched.
3. B1: Conversion-results table (source model → failure/check → fix →
   verification → residual limit → evidence path + machine) + proposed
   §6.3.1 paragraph, evidence-level tags on every row.
4. B2: plain-language Table 6.1 explainer + five-part-curriculum
   disambiguation table for Ben FIRST (§6a.3); then the per-task evidence
   audit with named gaps and the existing-results audit mapping each claim
   to its actual source/protocol. No table edits until Ben confirms.
5. B3: One connectome verdict: approved edge tables, programmatic edge
   diff, published-protocol baseline, ONE deletion/stimulus experiment with
   tolerances — or a documented blocker.
6. B4: Database live-schema check + one real linked example (record IDs +
   evidence) or an honest gap statement.
7. B5: Editorial items with fact-checks; overfull dispositions; the
   library-comparison item ends as ASK-THE-ADVISOR, never fabricated;
   the ethical subsection is untouched per §6a.2.
8. Updated CHATGPT_REPORT.md (in the restored copy) + supervisor verdict
   log per batch, all under `D:\GitHub\_zcode_backups\diss_campaign_20261024\`
   until Ben moves them back to EB475WS4.
9. BC (curriculum, this machine, MuJoCo): staged runs per §6a.8 —
   stages 1/2 air stepping, 3a/3b balance with harness then weaned,
   stage 4 harness walking with a literature-style harness (BWS %,
   sagittal constraint); stage 5 only if actually achieved. Per-stage
   traces + honest stage map of prior records.
10. BD: Deng + Li rebuilt from AnimatLab sources into base-Simulink,
    verified against AnimatLab headless runs and the numpy reference.
11. BIMG: redrawn leg/cup demo figures replacing the dissertation's
    current 60-sns images.

End every session by writing a one-screen status + open blockers file in
the campaign folder so Ben can follow from mobile.
