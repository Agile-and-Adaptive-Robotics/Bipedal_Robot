# Narrative polish and layer-specific verification

Completed in the current ProofFinal sources. No ZIP was extracted over working edits. No commits or pushes.

## Changes

- Chapter 6 now reports separate RG, PF, and MN interventions, with rhythm-generator and motor phase assessed independently. The two-cycle controls suppress both expected flexor burst onsets at every layer. RG stimulation shifts subsequent RG and motor timing by 0.211 cycles; PF stimulation and direct MN inhibition preserve the baseline schedule. These results concern the reduced Deng-derived circuit, not the full Gait2392 controller or the original rat mechanics.
- The full short-pulse and weak-current results are retained. Short pulses can truncate activity without deleting an entire burst. A one-second PF pulse shifts motor timing by 0.208 cycles while preserving the RG clock. It is explicitly distinguished from resetting the rhythm generator. Table D1 and Figure D4 show the protocols and representative two-cycle controls.
- Figure 6.7 overlays the fresh saved human-model trajectory with the Gait2392 inverse-kinematics reference, including stance/swing comparisons from foot loading. Both legs have ten retained cycles. No angle offset, amplitude scaling, mean subtraction, or fitted phase adjustment is applied. Each leg uses its own loading onset; this does not establish bilateral phase agreement. Table D3 gives absolute-angle RMSE: hip 28.03/24.65 degrees, knee 27.44/28.45 degrees, and ankle 8.71/7.38 degrees (right/left). These are newly calculated absolute-angle errors, not replacements for the existing mean-removed tuning metric.
- Narrative edits cover chapters 10, 15, 50, 60, 70, 80, and Appendix D. They clarify screening and citation terminology, define controller variants and gains, group the knee startup and loaded-arm explanations, place interface checks before full controller comparisons, distinguish completed reduced-circuit tests from planned full-network tests, and remove duplicated concluding prose. Chapter 7's sensory-control heading and proposed extensions now use direct language. Before/after sources and diffs preserve all revisions.
- The published muscle chapter, torque chapter, Appendix C, and torque-paper Main.tex are unchanged byte-for-byte during this pass. No biomimetic result, coefficient, torque curve, or measured point was changed. SADb counts remain 1,393 total app records and 943 entries in the Papers category.

## Fresh remote checks

Runs used the existing models and saved controllers, without optimization or retuning. Evidence is preserved under `Dissertation/Archive/flow_and_validation_20261028`; remote originals remain in `Code/MuJoCo_SNS/spinal/reports_20261028_publication_flow` on the heavy machine.

| Check | Result and scope |
|---|---|
| W2L fixed-controller feedback comparison | Reproduced intact phase 0.498 and displacement 2.239 m; feedback removal phase 0.223 and displacement 1.996 m. The removed-feedback run fails the specified supported-stepping phase criterion. |
| Human-model saved controller | Fresh score -145.0728022062, ten cycles per leg; model stance fractions 0.710/0.649 versus 0.610 reference. Joint overlay shows substantial hip and knee differences. |
| Deng-derived neural circuit | Sustained alternation, period 1.938 s; full 24-second accelerated/reference baseline voltage difference 0.0 mV. This verifies the translated implementation, not original-model quantitative equivalence. |
| RG/PF/MN stimuli | Initial weak sweep, six short strong controls, and three two-cycle strong controls retained. Each two-cycle run matches baseline neural and mechanical trajectories exactly before stimulation; baseline two motor onsets versus zero during intervention. |
| Rybak reduced circuit | All 13 existing checks passed in a fresh process. The first 240-second attempt timed out; the extended run completed in 302.5 seconds. No timeout is represented as a successful check. |
| Li MuJoCo translation | Incomplete replication: 0.040 m/s, eleven right-foot stance episodes, zero left-foot episodes. Appendix D distinguishes this failure from the earlier original-simulator result of 0.637 m/s. |
| Di Russo equations | Corrected Python/Numba equations agree to 3.33e-16 over 2,000 steps. This is internal implementation agreement, not equivalence to the author's native controller. |
| Di Russo original mechanical model | OpenSim model loads and initializes (nine active coordinates, eighteen Thelen muscles). Native custom-controller runtime remains unavailable, so published walking is not replicated. The rejected tumbling adaptation was not revived. |

The Deng-derived circuit follows the SNS-Toolbox demonstration described by Nourse et al. (2023), drawing on Deng et al. (2019). This pass does not claim replication of the separate 2023 Deng dissertation's full neuromechanical experiments.

## Verification and delivery

- Fresh five-pass local build: all return codes zero; no errors or undefined references; 314 pages. One existing 0.669-point overfull line remains in protected torque/Appendix C material. The newly introduced filename overflow was removed without changing results.
- New figures, tables, chapter transitions, and conclusion rendered and inspected. The neural stimulus plot was simplified to two columns to keep labels readable. A float barrier prevents the new controls from interrupting the archived-model paragraph.
- Live Overleaf: 315 pages, **Errors 0 / Warnings 0** after the final compile. The seven edited sources were read back from the live editor and compared with local sources; all match, with line endings normalized. Figures are in the established chapter folders.
- Root `Bolen_Dissertation.pdf` and `Bolen_Dissertation.zip` updated. Packaged source/support bytes verified against ProofFinal; 87 active graphics included. Source ZIP integrity check passed. The live-source ZIP download timed out, so verification instead used editor readback; no live ZIP was extracted.
- Paper source and delivery remain the previously verified 32-page main paper and 10-page supplement; this pass made no paper edits.

## Directive evidence

| Directive | Evidence in working source |
|---|---|
| Test RG, PF, and MN; resetting versus non-resetting | `60-sns.tex`, Layer-Specific Stimulation, lines 180–184; `95-AppendixD.tex`, RG, PF, and MN Stimulus Tests, lines 62–98. Exact raw values and complete cases: `Figures/95-AppendixD/Components/layer_stimulus_validation/stimulus_results_all.csv`; raw NPZ files and pre-stimulus equality checks in Archive. |
| Kinematics overlay against IK and stance/swing | `60-sns.tex`, Figure 6.7; `95-AppendixD.tex`, Direct Kinematics Overlay and Table D3. Curve CSV, generator, and methods in `Figures/60-sns/Components/gait2392_ik_overlay`. |
| Polish flow and jargon | Seven chapter diffs and `sources_before` in Archive; Chapter 6 interface/results order and Chapter 7-to-8 transition inspected. |
| Fresh relevant replication checks | Per-run logs, `rerun_summary.json`, extended Rybak status, source hashes, and scripts in Archive. Limits of Li and Di Russo reproduction are explicit. |
| Preserve protected/newer work | `protected_hashes.json` and `protection_verification.json`; no restoration from old downloads. |
| Sources and figures ready for Overleaf | `delivery_verification.json` and `live_source_verification.json`; live compile screenshot and figure screenshot in Archive. |

Build products are archived outside ProofFinal after verification. Figure ingredients remain in Components; active PDFs remain at the owning chapter's figure root. The remaining full-network stimulation studies and native Di Russo replay are research limitations, not completed replications.
