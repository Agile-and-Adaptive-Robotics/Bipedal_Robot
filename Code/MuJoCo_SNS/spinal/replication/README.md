# Replication program (Ben 2026-09-21): replicate Shinohara, Shevtsova,
# and Di Russo connectomes and results EXACTLY.

Status: FIRST PASS by ZCode — every file here is a draft that Ben will
edit. Nothing is authoritative until Ben signs it off.

## Why

The s3c-s3i verdict showed hand-built wiring does not produce bilateral
stepping. Replicating published, WORKING connectomes exactly removes
our wiring choices as a variable: if a published connectome walks in
our plant, we have a validated substrate; if it does not, the
difference is in the plant/rig, which is then a clean, publishable
finding.

## Targets (in build order)

1. **Di Russo 2023 (JNE 20 066006)** — MOST CONCRETE, first.
   Source: full text extracted (spinal\_dirusso_2023.pdf; rules in
   LIT_CIRCUIT_AUDIT.md section 10). Draft spec:
   `dirusso_connectome_draft.json`.
   - 9 muscles/leg, antagonist pairs per their Table 1.
   - CPG: two coupled phase oscillators (eq 6-7, heel-strike phase
     reset) + FIVE raised-cosine primitives (mu 0.1..0.9, sigma 0.2,
     eq 8) with SIGNED weights w_m,k (eq 9).
   - Reflexes: Ia mono+reciprocal via IaIN; II via two collateral INs
     (exc ag / inh ant); Ib autogenic via IbIN + IbIN mutual inh;
     Ib+ extensor force reversal; RC set (own MN + IaIN inh,
     antagonist-RC mutual).
   - NO state machine; receptor forms eq 11 (rIa = 65/200 *
     sqrt(max(0, v~)); rII = l~; rIb = f~; rf = foot contact force).
   - Foot: 3 spheres/leg (calcaneus r=5cm + two r=2.5cm at MTP line).
   - SUCCESS CRITERION (their result): walking 0.3-1.9 m/s; at
     1.17 m/s: step length 0.79 m, step duration 0.67 s, stance 60%,
     knee flexion peaks matching Bovi 2011 shaded ranges.
2. **Shevtsova 2026 (eLife RP107480)** — laminated rhythm generator.
   Source: full text held locally (see LIT audit). Draft skeleton:
   `shevtsova_laminar_draft.md` (V-class IN-mediated lamination;
   V1/V2b-mediated mutual inhibition, V2a RG->PF excitation, V0/V3
   commissurals). SUCCESS: half-center rhythm with laminated
   inhibition survives deafferentation at matching drive; period and
   duty in their reported range.
3. **Shinohara 2025** — interlimb/load adaptation model. Source:
   locate exact item in Zotero (search "Shinohara"); known anchor:
   force feedback excites the extensor center (their sec 4.2).
   Draft: `shinohara_draft.md`. SUCCESS: split-belt/treadmill
   adaptation behavior they report.

## How the replications plug into our stack

- Each spec JSON feeds the connectome editor/loader path
  (CONNECTOME.md); plant-specific bits (muscle list, foot geometry)
  are adaptation notes, clearly separated from connectome-exact parts.
- The phase-reset machine (pm_*) and contact-event machinery are the
  runner-side equivalents of Di Russo eq 7 and stay available.
- Per-replication validation figures go to replication\<name>\.
