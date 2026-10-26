# MID-SESSION RULING — inject into the RUNNING easteregg2 session now (Ben, 2026-10-24)

```
STOP AND RE-SCOPE — Ben's ruling, 2026-10-24 (overrides anything in the
running prompt that implies otherwise, including Batch BW's class-level
conductance plan):

1. MN input conductances and afferent feedback strengths are PER-TARGET:
   per-muscle, per-path, individual values. Collapsing them into family
   classes that share one strength is a waste of his time and is
   OVERRIDDEN. The 2026-09 "effective-strength doctrine" (per-instance
   conductances = GUI prefill, type-level only) is also OVERRIDDEN.

2. Re-scope Batch BW (and any build/tuning touching these pathways):
   Optuna solves INDIVIDUAL conductances per synapse — not class-level
   scalars. Seed literature values where they exist (Shinohara Table B.3
   per-MN PF->MN weights; Di Russo signed per-muscle W; eq 10/11 gains),
   Optuna the rest as individual bounded dims. Basin gate, sensitivity
   probes, and remove-and-retest are VALIDATION of the per-conductance
   winner, NOT justification to shrink the search space. Signs and
   topology stay fixed.

3. Inventory now: list every model/template in your scope whose MN inputs
   or afferent feedback share family-level strengths (base builders'
   G-family scalars, harvested templates whose mines dropped per-edge
   gains — "harvest cannot recover gain keys"). FLAG them for per-target
   rebuild; do not delete anything.

4. Mines/templates from now on must record per-instance values with the
   correct KIND: per-synapse SynAmp (uS) for chemical-synapse edges;
   per-path amplitude (nA) for afferent current-injection edges; MN to
   muscle couplings are activation mappings, not synapses.

5. The syn6 PF->MN pathway is already per-muscle (W-matrix Eq-18 mapping)
   — it is the pattern to follow, not the exception to remove.

Report the inventory + the re-scoped BW plan before continuing tuning.
```
