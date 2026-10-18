# W2L ROBOT WINNER → ANIMATLAB REALIZATION — stage report (2026-10-16, egg)

Overnight-campaign stage: **"AnimatLab realizer"** (phase 9, robot track). All
work ran ON easteregg2 (repo `D:\GitHub\Bipedal_Robot`, branch
KneeTestSetup_BenBo_stw, HEAD b497b67b) via ssh. NO commits, NO pushes; **no
file outside `reports_20261016\w2l_animatlab\` was touched** (proven:
`git status --porcelain Neuromechanical_Models` = empty post-run).

The ask: *apply the same L/R alpha+gamma drive + the highest-value feedback
fix from the campaign to the existing instrumented aproj/asim family in
`Neuromechanical_Models\`... generate a variant asim/aproj with the drive
inputs as stimulator objects on the L and R RG/MN populations and run it
headless (AnimatSimulator) to show the rhythm/stepping response. Label any
generated aproj 'headless-run only, GUI open-check pending Ben'.*

## 1. VERDICT

**The campaign winner (curr_s6_drive_robot_r060) transfers to the AnimatLab
Walker_2_Layer_CPG_BilateralRG model AS A SET, and the transfer is
self-consistent with the campaign's own mechanism story:**

| run (5.1 s headless, egg) | L/R RG flx bursts (−55 mV) | period (s) | antiphase (cycle) | note |
|---|---|---|---|---|
| **V0 control** (pristine copy) | 10/10 post-transient (10/11 raw incl. kickoff) | 0.467 / 0.462 | 0.481 | **reproduces the recorded 09-16/09-24 EB baseline** (11/11, 0.466 s, lag 0.474) |
| **V1 drive only** (gamma te/tf + alpha 0) | 2/1 | — | — | **LATCHES** (ext mean −57.4/−57.1 vs −58.8/−58.6 control; flexor suppressed after transient) |
| **V2 WINNER SET** (drive + comm + contact retune) | **13/13** (ext too: 13/13) | **0.361 / 0.361** | **0.4975** (n=12) | full-amplitude symmetric bursting, 23% faster rhythm, near-perfect antiphase; **deterministic 2/2 byte-identical charts** |
| **V3 feedback only** (no tonic) | 0 at −55 mV; **9/9 at −56..−58** | ≈9 cycles/4.7 s | — | rhythm PRESENT but SUB-THRESHOLD (peaks −55.6 mV) |

Read: gamma tonic ALONE kills the AnimatLab NonSpiking half-centers (the
documented 2026-09-15 tonic-latch doctrine holds for this neuron family);
feedback retune ALONE leaves the free-running rhythm just below burst
threshold; **together (the exact winner knob set) the model bursts at full
amplitude in clean antiphase** — the same "the tuned SET is what sustains the
weaned gait" result the MuJoCo campaign measured (wean report §3: defaults
die at 0.6 drive, no single pathway load-bearing). Contact-mediated stepping
continues in V2 (toe duty 21.9/38.8% vs control 15.9/38.7%; feet still
loading).

**Honest scope:** these are 5.1 s rhythm/stepping demonstrations at the
family's recorded-baseline scale — NOT the campaign's 30 s progression
criterion (that metric needs the MuJoCo yielding-walker protocol; this
AnimatLab model uses its fixed shelf harness, unmodified). No kine score is
claimed for the AnimatLab side.

## 2. BASE MODEL CHOICE (scoped to the night)

`Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\` — the bilateral-RG +
Shinohara-commissural + contact-adapter build (tools\SESSION_NOTES_20260916.md
state of record). Chosen because it is the **direct AnimatLab twin of the
MuJoCo winner net**: `spinal\w2l_mujoco\build_w2l_split_net.py` +
`build_w2l_aff_net.py` are the documented port of ITS wiring (M4 report;
2026-10-16 audit §1), so the winner's knobs map back onto named objects.
The Biped_2xCPG_wSubs walker was NOT given a variant tonight: the 2026-09-15
notes (`Biped_2xCPG_wSubs\tools\RG_contact_drive_notes_20260915.md`) already
proved its RG latches under ANY sustained tonic drive (matrix over
{2..60 nA} × G {5,0.5} × threshold {50,−55}), so applying the winner's
sustained te/tf there was known-dead; its contact-driven rebuild remains the
scheduled GUI-session work. The branch reference
(origin/AddingStepSensor_CoMorrow_stw Walker_2_Layer_CPG) is contact-driven
with NO tonic anywhere — its single-RG antiphase wiring is what the
BilateralRG build already subsumes bilaterally.

## 3. THE MAPPING (winner knob → AnimatLab edit)

Winner values: `w2l_wean\w2l_wean_best_params.json` (te=2.4, tf=3.0,
g_heel=0.58765, ctn_amp=0.71708, comm=0.41304; alpha=0; protocol knobs are
MuJoCo-test constructs with no AnimatLab counterpart).

| campaign knob | AnimatLab edit (asim + aproj in lockstep) |
|---|---|
| gamma te=2.4 nA | 2 tonic Stimuli "Gamma drive {L,R} RG ext (te=2.4nA)", 0→5.1 s, 2.4e-009 A |
| gamma tf=3.0 nA | 2 tonic Stimuli "Gamma drive {L,R} RG flx (tf=3.0nA)" |
| alpha_l/r = 0 | 4 Stimuli "Alpha drive {L,R} {Hip,Knee} MN ext (0nA)" — inert at 0 A, GUI-visible knobs for Ben |
| comm=0.41304 | "c1 Commissural Inhibit" SynAmp 2.749→**1.1354**; "V3 Commissural Excite" 0.1→**0.0413**; NEW type "RG Excite (comm drive, campaign)" SynAmp **1.1354** repointing the 4 RG→{c1,V3} input connexions (the shared "RG Excite" type — used by the RG half-centers themselves — untouched) |
| ctn_amp 0.717 of 2.0 | 4 ContactCount adapter Gain C **20 → 7.1708 nA** (design-intent × campaign ratio ×0.3585 — see unit-bug below) |
| g_heel 0.588 of 0.5 | NEW type "Afferent HC Excite (Contact, campaign)" SynAmp 0.01→**0.011753** repointing ONLY the 20 contact→extensor-layer connexions (the shared "Afferent HC Excite" stays on the 24 Ia/II/Ib afferent links) |
| g_ib=0.112 / ibnA=1.321 | **NOT transferred** — the campaign's own sensitivity probe at the winner shows the Ib pathway non-load-bearing there (removal Δobj = 0, wean report §3), AND the Ib afferents share the "Afferent HC Excite" type with the Ia/II flexor afferents that have no campaign analog |
| toe_m ×2.2198 | **N/A** — no dorsiflexion path exists in this AnimatLab net (same documented reduction family as the MuJoCo port) |
| cref=46.37 N | **N/A** — AnimatLab adapters are ContactCount-based, not force-based |

Mechanism check that the comm edit reached the circuit: in V2 the L/R c1 and
V3 cells quiet down (max −56.1/−55.9 mV, never crossing −55; control maxima
−53.7/−52.6) — the commissural path is measurably weaker, exactly the
campaign's "the two RGs need LESS mutual inhibition once heel feedback paces
them" (wean report §5).

**Structural gap documented (Ben's next rules pass):** the AnimatLab contact
wiring is ipsilateral-extensor excitation only; Ben's rules crossed
heel→contralateral RG-flx-IN reset exists only in the MuJoCo M5 build.

## 4. UNIT BUG FOUND (flag for Ben)

The 4 ContactCount adapters shipped with `<C>20</C>` in the asim — plain SI
**amps** per contact count — while every other adapter in the file uses
amp-scale values (e.g. 1.741e-009); the aproj mirror reads
`<C Value="20" Scale="nano" Actual="1.741e-009"/>` (Value and Actual
inconsistent; the Actual is stale from the cloned II-adapter template).
`build_contact.pl` evidently intended **20 nA**. As shipped, contact input
saturates the contact neurons (a binary stance switch) — which is how all
recorded ground runs behaved. Tonight: **V0 keeps the shipped C=20 (faithful
control); V2/V3 carry 7.1708 nA** = 20 nA design intent × the campaign's
encoder ratio. Even the corrected 20 nA is ~28× the MuJoCo winner's 0.72 nA
encoder scale — only RATIOS transfer across the two stacks (their contact
architectures differ). Ben should re-export the stock aproj's adapters from
the GUI (which recomputes Actual) or fix build_contact.pl when convenient.

## 5. FILES (all under `spinal\reports_20261016\w2l_animatlab\`)

- `V0\ Walker_2_Layer_CPG_BilateralRG_V0_control.asim` — byte-identical copy
  of Ben's ground asim (proven `fc /b`), + chart .txt byproducts + analysis.json
- `V1\..._V1_drive.asim` — + 8 stimulators (gamma/alpha) only
- `V2\..._V2_winner.asim` — stimulators + comm + contact retune (**the
  realization**) — deterministic 2/2
- `V2\..._V2_winner.aproj` — aproj mirror of V2 (same edits, GUI-visible
  objects; drawings/pages untouched by construction, per-page CDATA count
  check passed). **HEADLESS-RUN ONLY, GUI OPEN-CHECK PENDING BEN** — never
  saved from the GUI, not opened in AnimatLab2 tonight.
- `V3\..._V3_feedback.asim` — comm + contact retune, no tonic gamma
- `build_winner_variant_20261016.py` (deterministic GUID prefix d51a,
  XML + page-CDATA verification built in), `analyze_run_20261016.py`,
  `probe_crossings.py`, `edits_log.json` (machine-readable edit list)

V1/V3 are throwaway asim probes (A/B evidence); V2 is the only aproj variant,
per the no-asim-only-backdoor rule for model changes.

## 6. REPRODUCE (egg, cmd)

```
cd /d "D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin"
AnimatSimulator.exe "D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\w2l_animatlab\V2\Walker_2_Layer_CPG_BilateralRG_V2_winner.asim"
D:\Anaconda\envs\myo\python.exe D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\w2l_animatlab\analyze_run_20261016.py D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\w2l_animatlab\V2
:: rebuild everything from Ben's pristine model files:
D:\Anaconda\envs\myo\python.exe D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\w2l_animatlab\build_winner_variant_20261016.py --src D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG --out D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\w2l_animatlab
```

Analyzer definitions: chart volts→mV; trailing all-zero rows popped;
rhythm = flexor up-crossings of −55 mV (the 09-18 session rule); antiphase =
median L-onset phase within consecutive R onsets. The OSG viewer fails to
create a window over ssh (log lines retained in analysis folders' console
output) — the simulation itself runs to SimEndTime 5.1 and exits 0.

## 7. NOT RUN / NOT CLAIMED (honesty list)

- **No GUI open-check** of the V2 aproj (AnimatLab2.exe exists on egg but a
  window-station-less ssh session cannot produce a trustworthy pass/fail) —
  labeled pending Ben. Suggested check: open V2 aproj, expect zero Error
  dialogs, 10 stimuli in the tree (2 stock kickoffs + 8 new), then re-export
  the asim from the GUI if adopting.
- No 30 s / longer-horizon persistence run (kept the family's 5.1 s baseline
  scale for comparability); no COM-progression metric on the AnimatLab side;
  no basin gate (no AnimatLab harness exists for it — determinism 2/2 shown
  instead).
- Ib and toe chains not transferred (reasons in §3); no crossed
  heel→contralateral-InF edge added (topology change beyond tonight's scope).
- The winner's harness/protocol knobs (legcoll/sink/rig/walker_yield) have no
  AnimatLab counterpart — the shelf harness was left as shipped.
- Biped_2xCPG_wSubs: no variant (known tonic-latch; §2).

— AnimatLab realizer, 2026-10-16
