import os

ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022"

CUT2 = """

## Round-2 cuts (2026-10-22, after principal review: "selection, compression, and explanation — not another expansion")

60-sns (10,851 -> 6,794 narrative; captions 1,069 -> 1,096; table bodies 316 -> 463;
combined reading matter 12,236 -> 8,353, -31.7%):
- Four-tool-families paragraph reduced to a one-sentence pointer; tab:sns_capability
  is now the primary tool description (round 1 had table + full prose restatement).
- "The standard of proof..." one-liner deleted (duplicated the verification framing).
- AnimatLab-organization paragraph: instrumentation/bilateral-rebuild detail folded to
  three sentences; figure set unchanged.
- MuJoCo conversion: repair-script narration compressed; "Python reproduces MATLAB"
  parity sentence, 0.9 s rhythm-generator detail, and the Selionov tonic-drive aside
  CUT (logged; the parity fact remains implied by the audited-pipeline claim).
- Simulink library: block inventory to one sentence; CAD hinge-mate defect kept
  verbatim-precise (relevant limitation), URDF-trial detail cut.
- Ports, cross-engine methods, spiking methods: each ~40-55% shorter; protocol
  details retained in results or Appendix; sub-stepping/residual numbers moved out
  (Appendix tabulations remain).
- Results: Simulink demos + ports compressed ~45% (10 nA/20 ms pulse, 96%-recruited,
  block-count breakdown cut); AnimatLab validation per-model SDs/duty/threshold
  detail cut; gait-sweep sub-artifact sentence (-8816/-8761) cut; deletion-stimulus
  per-perturbation scores and pelvis-tilt numbers cut (table + prose summary kept).
- Spiking results: unit-level residuals moved to "tabulated residuals"; tab:spike_beer
  DELETED (six rows folded into one prose sentence: sag 2.02 vs 2.00, off 7.85,
  modulation within 0.02 deg); hybrid-subnetwork paragraph halved.
- Discussion: SIX subsections -> THREE. sns_bpa_integration + bpa_walker +
  sns_controller_comparison merged into "The BPA-Actuated Robot and the Walking
  Status" (label sec:sns_controller_comparison preserved; sec:bpa_walker and
  sec:sns_bpa_integration labels retired, no external references — verified);
  sns_tools_summary subsection DELETED (repeated the capability table; replaced by a
  two-sentence closing handoff). Walking status now stated once.
- CONTRADICTION FIXED: curriculum paragraph said the -25 fallback "sits below every
  genuine gait" — wrong direction; it sits numerically ABOVE every genuine gait
  score (genuine gaits -144 to -267; -25 > -144), which is exactly why it is a
  sentinel excluded from score comparisons, and the text now says so. Matches
  tab:walker_ablations/tab:gait_sweep captions.

50-database (2,588 -> 2,063 narrative):
- Opening restructured DATABASE-FIRST: the SADb is the subject of sentence 1; the
  controller-requirements motivation compressed from ~200 to ~130 words.
- Methods: import-batch/recovery counts cut (374-in-fifteen-batches, 47-of-55 PDF
  recovery, 24-of-70 abstract recovery -> one clause); prune three-stage enumeration
  merged into one sentence (95/290/4+47 stage counts cut, totals 436/382/54 and the
  351-protected clause kept); "each approved by me before execution" dropped
  (decision authority already stated once, polished, in the division-of-labor
  paragraph).
- CONTRADICTION FIXED: "DOI as the unique identifier of every record" vs 550-of-561
  with DOIs -> "unique identifier wherever a record carries one".

App-first ruling applied (Ben, 2026-10-22): the packaged explorer is the database's
primary form; Airtable is the editable working base; every explorer build embeds the
full 943-record audited corpus including the 382 papers pruned from the working base
(verified against export/sadb_export.json archived flags and app check_app.py:
"papers embedded: 943"). SADb_audit/app/README.md roadmap line updated to match.

bolen_design_2020 citation REVERTED from 60-sns (Ben: it is a poster schematic, not
a paper, and its connectome is not the current one). No new citation added there.
"""

MOVE2 = """

## Round-2 (2026-10-22) paragraph accounting — SNS

Deleted outright (content absorbed by tab:sns_capability or the merged discussion):
- "The standard of proof for each tool..." (1 sentence).
- tab:spike_beer (table; numbers preserved in one prose sentence).
- sec:sns_tools_summary subsection (2 paragraphs) — its unique surviving content is
  the closing handoff sentence pair ("What the tools are for is the program they
  define next..."; "Every tool is reusable...").
- The controller-comparison opening paragraph and BPA-integration subsection were
  MOVED (not deleted) into the merged "BPA-Actuated Robot and the Walking Status"
  subsection (sec:sns_controller_comparison), each compressed.
Labels retired: sec:sns_bpa_integration, sec:bpa_walker, sec:sns_tools_summary
(no external references; verified by grep across all chapters). All externally
referenced labels preserved: sec:prior_walker, sec:preliminary_sim_methods,
sec:spiking_methods, sec:sns_tuning_methods, sec:preliminary_sim_results,
sec:gait_library_results, sec:walker_variant_results, sec:disc_spiking,
sec:sns_controller_comparison, sec:sns_gaps, sec:sns_curriculum.
Every verified number removed from prose in round 2 is listed in CUT_MAP round-2
section; none were changed, and each remains traceable to its source table, figure,
or Appendix tabulation.
"""

for fname, text in (("CUT_MAP.md", CUT2), ("SNS_MOVE_MAP.md", MOVE2)):
    with open(os.path.join(ROOT, fname), "a", encoding="utf-8") as f:
        f.write(text)
    print("appended round-2 section to", fname)
