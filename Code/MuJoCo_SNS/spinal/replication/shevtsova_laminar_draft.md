# Shevtsova replication draft — UPDATED 2026-09-22 (full text in hand)

Source: Shevtsova et al. 2026, eLife 107480, "Linking spinal circuit
reorganization to recovery after thoracic spinal cord injury"
(DOI 10.7554/eLife.107480). Local full text:
`spinal\shevtsova_2026_fulltext.txt` (123,725 chars, from Ben's local
Zotero storage 6LM33WX6 .zotero-ft-cache; the PDF itself is not synced
to the servers — cite the local cache).

## What the model IS (from the full text)

The RAT adaptation of the laminar quadruped locomotor circuit model
lineage: Danner et al. 2017 (eLife 6:e31050) and Zhang et al. 2022,
adapted to rat, then subjected to two thoracic SCI conditions
(lateral hemisection; contusion) to link circuit reorganization to
the observed gait changes. This is the SAME model family our RG/PF +
V-class architecture copies.

Population inventory (from full-text keyword evidence):
- Bilateral RG-E / RG-F half-centers WITH persistent NaP (12 NaP
  mentions; our HC class matches).
- Laminated mutual inhibition via the V1/V2b-mediated INs (our
  InE/InF correspond — verify naming against their figures).
- Commissural pathways: V0 (51 mentions — V0D/V0V splits), V2a
  (7), V3 (15) — the left-right coordination set, richer than our
  coarse c1/V3 pair.
- Speed-dependent gait expression (walk/trot transitions) from
  descending drive — the behavior our DRIVE knob approximates.
- 28 RG-F / 4 RG-E mention hits: the text focuses on flexor-side
  reorganization after SCI.

## Replication plan (Ben edits before build)

1. EXACT SOURCE: the rat model files. eLife publishes model code —
   pull the Danner/Zhang-lineage rat model (ModelDB / GitHub link in
   the paper's data availability; TO EXTRACT from the fulltext
   "Data availability" section — next pass).
2. Populations: map 1:1 to our SNS vocabulary (NaP HCs, V-class INs,
   PF layer, MN pools; quadruped -> biped adaptation of the
   forelimb/limb-1 populations is the ONE structural decision Ben
   makes).
3. Connectome: their figure tables -> edge list via the block editor
   (load the bilateral literature template as the starting point, then
   Ben edits to their figure).
4. SUCCESS CRITERIA (fill after 1):
   - [ ] intact model: speed-dependent stepping, walk-trot expression
   - [ ] hemisection: their reported gait asymmetries reproduce
   - [ ] contusion: their reported deficits reproduce

## Standing notes

- The V0D/V0V split (inhibitory/excitatory crossed inhibitory paths)
  is the structural piece our current c1/V3 lacks — the s3g matrix
  (crossed-cut freed the frozen leg) makes this replication directly
  diagnostic for OUR architecture failure.
- NaP half-center parameters: we already carry the fixed-tau-h class
  with the Deng/Di Russo parameterization; Shevtsova's parameter set
  is the literature check on ours.
