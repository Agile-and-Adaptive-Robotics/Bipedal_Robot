# Replication program (Ben 2026-09-21): replicate Shinohara, Shevtsova,
# and Di Russo connectomes and results EXACTLY.

Status: FIRST PASS by ZCode — every file here is a draft that Ben will
edit. Nothing is authoritative until Ben signs it off.

## Why

The s3c-s3i verdict showed hand-built wiring does not produce bilateral
stepping. Replicating published, WORKING connectomes exactly removes
our wiring choices as a variable: if a published neuromechanical model
(connectome + musculoskeletal model, coupled) walks here, we have a
validated substrate; if it does not, the difference is isolated to our
musculoskeletal model/rig, which is then a clean, publishable finding.
NOTE (Ben): these are NEUROMECHANICAL MODELS, not "plants" - the
musculoskeletal side is part of the coupled system its reflexes loop
through, not passive plumbing. Replication therefore carries the
paper's musculoskeletal parameters (Di Russo: Delp skeleton, 1.8 m,
75.16 kg, Thelen muscles - gait2392 lineage, which our conversion
shares) AND its coupled controller-optimization loop, not a bolt-on.

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
   Source: full text held locally (see LIT audit). Draft:
   `shevtsova_laminar_draft.md` (FULL EXTRACTION 2026-09-23: Table 1
   weights verified, lumbar-only biped reduction) + machine-readable
   rules `shevtsova_rules.json`. SUCCESS: half-center rhythm with laminated
   inhibition survives deafferentation at matching drive; period and
   duty in their reported range.
3. **Shinohara 2025** — interlimb/load adaptation model. RESOLVED
   2026-09-23: bioRxiv 2025.11.11.687930 (cat hindlimb + Danner-lineage
   CPG + per-muscle afferents INTO the centers), full text local
   (`shinohara_2025_biophiv687930_fulltext.txt`). Draft:
   `shinohara_draft.md` (full extraction incl. Tables A.1/A.2/B.3 +
   eqs 9-11 afferent wiring) + `shinohara_rules.json`. SUCCESS:
   split-belt/treadmill adaptation behavior they report.
4. **Rybak 2006a/2015/2024** — two-level RG/PF rhythm model. Draft
   `rybak_draft.md` + `rybak_rules.json` (2026-09-23; 2006a Table 2
   extracted from PMC1890439 — full text now local as
   `spinal\lit_rybak2006_fulltext.txt`; companion 2006b is paywalled,
   NOT read).

## How the replications plug into our stack

- Each spec JSON feeds the connectome editor/loader path
  (CONNECTOME.md); plant-specific bits (muscle list, foot geometry)
  are adaptation notes, clearly separated from connectome-exact parts.
- The phase-reset machine (pm_*) and contact-event machinery are the
  runner-side equivalents of Di Russo eq 7 and stay available.
- Per-replication validation figures go to replication\<name>\.
