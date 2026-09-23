# CONNECTOME.md — Ben owns the wiring, the machine tunes

Established 2026-09-21 after the s3c–s3i verdict: output-level tuning
cannot fix a mis-wired circuit, and architecture decisions are Ben's.

## Workflow

1. **Open `connectome_editor.html`** (double-click; runs in any
   browser, no server).
2. Each rule row = one literature connectome rule (Deng A6 / Di Russo
   rules 1–5):
   - **enabled** — rule present in the circuit or removed entirely.
   - **hops** — how many interneurons sit between the afferent and its
     target. 0 = DIRECT edge (legacy/short-cut wiring — the thing that
     broke the rhythm); 1+ = the literature interneuron layer.
   - **gain** — connection strength (sign is fixed by the rule's
     biology; shown next to the field).
3. **Export** writes `connectome_gains.json` beside `build_network.py`.
4. Any runner run (or study) applies it BEFORE building the network:
   enabled rules set their gain keys, disabled rules zero them, and any
   rule with hops >= 1 turns on `full_rules` (the interneuron-layer
   topology: IIX/IIIN/IBIN populations, IaIN/IbIN mutual inhibitions,
   heel/toe via InE/InF lamination).

## Rules currently in the editor

| id | rule | default hops | gain key |
|----|------|--------------|----------|
| ia_homo | Ia -> agonist MN (monosynaptic) | 0 | ia_to_mn |
| ia_recip | Ia -> IaIN -> antagonist MN (+ IaIN-IaIN mutual) | 1 | ia_to_antagonist |
| ii_exc | II -> IIX -> agonist MN | 1 | ii_to_mn |
| ii_inh | II -> IIIN -> antagonist MN | 1 | ia_to_antagonist |
| ib_auto | Ib -> IBIN -> same MN (+ IbIN mutual) | 1 | ib_to_mn_inh |
| ib_rev | Ib+ stance-gated extensor reversal | 1 | ib_group_exc |
| heel | heel SN -> IN -> RG via InE/InF | 1 | heel_rge |
| toe | toe SN -> IN -> RG via InE/InF | 1 | toe_rge |
| ib_load | stance-Ib group IN (LBIN) -> RG-E | 1 | ib_rge |
| rc | Renshaw set | 1 | renshaw |
| cross | crossed commissurals c1/V3/contra | 1 | c1_gain |

## Division of labor

- **Ben**: every wiring decision (rules on/off, hops, gains). The
  diagram is generated from these decisions, so the figure can never
  disagree with the circuit again.
- **The machine (ZCode)**: implements the rule semantics faithfully,
  runs the sims/studies, tunes numeric parameters ONLY, and reports
  honestly.

## Notes

- **Block editor**: `connectome_block_editor.html` (v2, Ben-requested
  block coding) — palette of prebuilt types (V/C INs, SNs, HC RG/PF,
  reflex INs, Renshaw, MN, muscle), drag-place, SHIFT-click pairs to
  draw a synapse, per-synapse sign/gain/tuning-tag, Save/Load JSON,
  "Load Deng A6 motif" seed. `w2l_equivalent_draft.json` loads the
  Walker-2-Layer BilateralRG equivalent (contact-driven RG per the
  09-15 notes) for editing. Delete/Backspace deletes the selection;
  undo/redo + arrowhead-on-body are queued (v2.1).
- `hops = 0` on any rule reproduces the legacy direct wiring; the
  09-21 deafferentation matrix
  (`_deaff_matrix.py` -> `_deaff_matrix_out.json`) is the evidence
  pack for why the interneuron layer + crossed-communication care
  matter: cutting the crossed paths alone let the frozen leg step
  (7 cycles) for the first time.
- The builder raises on unknown rule ids (editor/builder version
  mismatch = loud failure, not silent rewiring).
