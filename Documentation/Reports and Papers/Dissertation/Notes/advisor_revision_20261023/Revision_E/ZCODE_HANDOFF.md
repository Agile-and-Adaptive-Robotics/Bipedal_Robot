# ZCode handoff — Revision E

Read `Revision_E/REPORT.md`, `review/parameter_plumbing_audit.json`, `review/self_scoring_matrix.json`, and Ben's replies before changing prose. The chapter's point is to establish what the tools let us test about neural control, feedback, and mechanical differences. Completing a diagram, import, search, or build is not acceptance of that scientific claim.

## First: recover the exact running source

Run `hostname` first and use only that machine's AGENTS.md entry. Paths below are repository-relative. On easteregg2 the checkout is normally `D:\GitHub\Bipedal_Robot`; do not copy desktop environment paths. Preserve local/uncommitted work and record source hashes, backend versions, model hash, actual command, and output hashes. The fetched desktop baseline is a988508a. A newer easteregg2 runner may explain the discrepancy; do not overwrite it with this baseline.

Read the `dissertation`, `animatlab`, `mujoco`, `matlab`, `neuroscience`, `sns-toolbox`, and `opensim` skills from the local ZCode_Skills repo/junctions. Use relevant configured plugins where they provide direct source or model access. Read `spinal/DESIGN.md` before runner/builder edits. Preserve the s3k anchor and its existing contracts; separate W2L and Di Russo paths. Do not equate an adaptation with original-paper reproduction.

## Priority 1: parameter plumbing before optimization

At the reviewed commit, `reports_20261016/w2l_wean/_w2l_wean_eval.py:evaluate` keeps only keys already present in `T.KV`. The selected w25 file has 61 fields; **55 are dropped**, leaving six legacy knobs. `w2l_mujoco/test_w2l_ground.py:run_phase` constructs family-level gains and does not forward the individual per-target fields, despite their existence in `w2l_cpg/build_w2l_net.py`. The sensitivity file's many exact-zero effects therefore cannot be treated as evidence that those pathways are unimportant.

1. Confirm the exact source used on easteregg2. If its loading differs, save that version with its original outputs first.
2. Make each intended key explicit in accepted input, builder forwarding, saved winner, and reload. Reject unknown search keys loudly; never silently discard them.
3. For every searched key, perturb that key and prove the expected compiled edge conductance changes, while unrelated edges retain their values. Record edge IDs, signs, units, and source-specific conductance interpretation. Distinguish encoding gains from synaptic conductances. Shared L/R parameters must be labeled reductions.
4. Save the effective compiled configuration, not only the intended JSON. Reload in a fresh process and compare signatures. Check default settings against the retained baseline.
5. Replay the existing w25 and w10 candidates with the exact old effective configuration and corrected intended configuration. Do not silently replace earlier records or claim earlier scores apply to corrected dynamics. Only then decide whether additional search is needed.

## Priority 2: deletion and perturbation evidence

`reports_20261025_persyn/_poster_panels.py:cycles` always uses `heelF[:,0]`; the right count calls the same function on the same matrix. Give each leg its own trace. Count after the settled window with the same threshold/merge conventions as the evaluator. Its `kneeL_deg` uses hip column 0 and converts degrees again; log q is already degrees in JNTS order, knee_L index 1. Correct names and units rather than interpreting the large numbers physically.

Before every deletion, verify the target edges exist and that their effective compiled weights become zero. Restore the baseline between cases. Several current deletion keys do not reach the builder. Do not infer redundancy from an unapplied deletion. Ib's `ibnA` encoding removal is distinct from deleting specified synaptic edges.

Recover saved raw logs first. If absent, rerun only the necessary baseline/deletion/stimulus/push conditions after the plumbing checks. Record bilateral contacts, joints, support force, forward displacement/speed, phase, duty, fall/finite status, and post-perturbation windows. Define recovery tolerance and horizon before interpreting the push. State the 0.15 s pulse duration, force direction, stimulus amplitude, and 1 s duration. Show the actual rather than nominal harness contribution. Do not treat nominal 10% wsup as measured 10% total support.

The mass-change evaluator independently counts both legs and its ten records remain useful, but match them to the effective controller/source. Its driver inherits other runner defaults; audit `pushN` and every perturbation switch, including whether a default 40 N pulse was present. Record unperturbed controls explicitly. Preserve body/inertia changes and required physiological joint/contact gates. This is W2L, not the 92-actuator BPA robot model.

## Priority 3: explain what the additional work establishes

- Deng-derived Simulink: ten bursts / 1.938 s period after one pulse is an intrinsic-rhythm test under fixed tau_h. The log also shows about 63 mV max waveform differences. Diagnose alignment, solver, and equations before calling it a dynamic replica. Do not replace source tau_h(V) in other literature models without labeling that change.
- Full 410-neuron comparison: all final voltages at 0.3 s and nine short traces pass; extend to cycle statistics and closed-loop behavior separately. Short agreement does not prove long-term transfer.
- Li: 58/98 topology with graded adaptation replaces original spiking stance cells. Scripted antiphase contacts and correlation -0.271 are an adaptation check; verify the original closed-loop experiment separately.
- Air-step: fixed pelvis, no contact, eight BPA elements; mechanics about 0.31 s versus RG 1.938 s. Establish neural-to-mechanical entrainment and feedback effects before calling a curriculum stage or gait validated.
- Di Russo: committed stage-2 search has 120 failed trials and best objective -100. Preserve that negative result. Inspect validity gates, muscle mappings, balance and feedback equations, receptor delays, stance/contact, and stage-1 excitation fit before spending more optimization. Do not weaken gates to manufacture a pass.
- Hunt/Fletcher library comparison remains pending identity/source from Ben. Compare editable neuron/synapse coverage, sensory/actuator coupling, mechanical interfaces, and actual verification. Avoid an unsupported novelty claim.
- Newly produced local v4 demo figures are not in E. Preserve their files and Ben's figure instructions. Knee-reflex v4 reports max 3.6° discrepancy; label approximation and show timing/solver evidence. Do not call those twins exact.

## Work allocation and acceptance

Ben prefers GLM-5.3 / GLM-5.3-Flash for code, modeling, testing, and optimization, and ChatGPT Astra for final prose. Use the smallest evidence-recovery and plumbing checks first; heavy runs belong on easteregg2. This handoff authorizes no new external message, Overleaf edit, commit, or push. Ben will initiate ZCode on the other machine.

Return one evidence packet mapping each claim to source version, model, neural variant, effective parameter file, protocol, raw outputs, criterion, result, and limitation. Include failed runs and substitutions. Then update the local E draft with precise results. Keep repeated clarity concerns unaccepted until Ben/advisor can explain the goal and result back without build-order jargon. Every dissertation revision ends in a clean build plus line-based compliance evidence.
