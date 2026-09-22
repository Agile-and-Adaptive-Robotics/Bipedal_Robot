# Shevtsova replication draft (FIRST PASS — needs the full-text pass)

Source: Shevtsova et al., eLife RP107480 (2026), "Laminar organization
of the mammalian locomotor network" — full text held locally per
LIT_CIRCUIT_AUDIT.md; ALSO the earlier Shevtsova/Rybak laminar RG
papers. THIS FILE IS A SKELETON for Ben to edit; the full-text
extraction pass has not run yet.

## What to replicate

The laminated mammalian rhythm generator: rhythm and pattern formation
with ALL mutual inhibition routed through V-class interneurons, at
published conductances, reproducing (a) alternation with correct
phase relationships, (b) their reported sensitivity analyses.

## Population skeleton (per side)

- RG-E, RG-F: persistent-Na half-centers (we have the class).
- V1 + V2b: inhibitory INs mediating RG-E->RG-F and RG-F->RG-E
  (our InE/InF correspond; verify V-class assignment and conductance
  split against the paper's figures).
- V2a: excitatory RG->PF and RG->MN communicaation (our PF excitatory
  drive corresponds; check gain structure).
- V0 (V0D/V0V) + V3: left-right commissural pathways (excitatory V3;
  inhibitory V0D; our c1/V3 are coarse stand-ins — the paper gives
  the proper split).
- IaIN, IbIN, II INs, RC: reflex layer as in Di Russo rules but with
  Shevtsova's conductances.

## Success criteria (fill from paper after full-text pass)

- [ ] period at stated drive
- [ ] duty and phase relationships (E/F, left/right)
- [ ] deafferented rhythm persistence
- [ ] their key ablation figure reproduced

## TODO before build

1. Full-text extraction of the held PDF (same pipeline as Di Russo:
   pypdf + section sweep).
2. Pull their figure tables of conductances (they publish them).
3. Map to our SNS population vocabulary; flag any population our
   toolbox cannot express.
