# ZCode handoff: evidence and integration work after the local prose pass

Read CHATGPT_REPORT.md first. Ben authorized local editing and asked for this handoff; Overleaf remains read-only. Do not upload or replace chapters from this folder without a new live-source comparison. Do not commit or push. This working copy is Revision_D; its changed chapters are Introduction, Background, Database, and SNS. Use review/changes.json and changes.diff to merge targeted edits, not whole files. Frozen before copies are in review/before.

## Routing and acceptance

Ben wants GLM-5.3 / GLM-5.3-Flash for testing, modeling, and optimization through his ZCode setup. Use the existing configured access, not a newly created provider integration. Final prose polish is for ChatGPT Astra, limited to evidence-checked changed sections. Do not spend a premium-model pass on broad repository rediscovery or loop indefinitely polishing the same paragraph. Do not optimize before verifying transcription and baseline behavior.

## Priority 1: provide actual findings, not more abstract prose

1. **Conversion checks (SNS §6.3.1).** Inspect repair_robotbody_muscles.py and the generated/source model comparisons. List the actual muscle-path/joint issues corrected, if any, and the exact checks that passed. Return source file/function, before/after model property, a rerunnable check, output, and proposed two-sentence result. The advisor's wrong-origin/left-hand-axis examples were hypothetical, not findings. Do not label a procedure successful from its existence alone.
2. **Five-system status (SNS §6.2 / Table 6.1).** For each of the five tasks, give implemented model, neural controller, mechanical model, actuation, input protocol, output evidence, and explicit remaining gap. Confirm bilateral PF organization and contributors for the old 'for each leg?'/'We?' comments. The local pass retained the captured source's I/Morrow attribution; it did not independently establish authorship history. The fifth task is not complete merely because the SolidWorks mechanism imports. Check the rigid ankle limitation and whether neural/BPA/mechanical coupling actually ran.
3. **Existing results first.** Audit CROSS_PLATFORM_REPORT.md and its referenced logs, BRIDGE_REPORT.md, E1 reflex and E2 open-loop experiment scripts, and W2L matched-circuit records. Separate open-loop 92-actuator drive, closed-loop reflex exchange, driven rhythm-core agreement, and full neuromechanical replication. Report which actual source/protocol matches each claimed result. Do not rerun tests already adequately evidenced unless a specific mismatch requires it.

## Priority 2: one exact connectome and one published protocol

Read Code/MuJoCo_SNS/spinal/DESIGN.md before touching runner.py/build_network.py/params.py. Preserve the default census and bit-exact anchor contract; use fresh processes for relevant verification. Keep MuJoCo 2.3.7 and the recorded step through this work. Heavy runs belong on easteregg2; Simscape runtime work belongs on the laptop.

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
- Background's remaining ethical-case and scientific-validation claims need a scoped fact check: broad replacement of animal experiments, unsupported validation-chain claims, and claims that all interventions become routine should not be polished into stronger assertions. The present pass fixed repeatability and throughput wording only.
- Background Haeufle robustness wording: check the actual paper before replacing a technical result with a plausible paraphrase.
- Earlier-library comparison: identify the actual Hunt/Fletcher libraries meant by the advisor before comparing them. Do not substitute a different library or remove a valid requested comparison because the citation is missing.
- Build support outside the five captured chapters came from Revision_C. Reconcile the unedited chapters, main.tex, bibliography, and figures with a fresh complete Overleaf export before any upload. Bibliography was copied byte-for-byte, not re-encoded. The existing Flx20mm_meetsHuman.pdf had to be copied from tmp/overleaf_live_20261023 into this working copy for compilation.
- The build has 26 overfull boxes, predominantly outside this targeted prose work; locate and disposition them. Do not claim clean layout from zero compilation errors.

## Return contract

For every handled matrix row, return: original request, section/line and exact before/after, evidence and command/log, result, remaining limit, and self-score. A clarity revision can earn 3/4 after these checks; 4/4 requires Ben/advisor acceptance. Repeated concern reopens the row. No GO based solely on a build or a populated checklist. Update the matrix and CHATGPT_REPORT.md, preserving old scores and their evidence. No administrative report text belongs in the dissertation.
