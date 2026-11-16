# Independent historical W2L review

Read-only reanalysis of retained primary recordings and scoring definitions. No simulations, fitting, plotting, or canonical-source edits were performed. Reproduce using `D:\Anaconda\python.exe HISTORICAL_W2L_audit.py`. Full path/hash bindings and computed results are in `HISTORICAL_W2L_RESULTS.json`. The associated historical human-model IK review is `HISTORICAL_IK_REVIEW.md`, with a separate independent raw-array recomputation.

## Primary sources and identities

Repository prefix: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\`.

| File relative to prefix | SHA256 |
|---|---|
| `reports_20261016\w2l_wean\_w2l_wean_eval.py` | `bad355e83c31063b393af79465dea21e9ef5a3f0b91b9c9071b0a4b0c795fd76` |
| `reports_20261027_codex\w2l_verified_config.json` | `868cb653bea2bd91ef8c088fb52562f53520233d187c981f7fe7f15a35cb4787` |
| `reports_20261027_codex\trial458_d30.npz` | `53dd98df62d00d920e8ea4e329f7e4b9f7558d36bc0d155459c133d07fb45875` |
| `reports_20261027_codex\trial458_d30_ws0.npz` | `b7c316143b78847bfd6ae42af0b230039366844fe9703172e870703a07862b12` |
| `reports_20261027_codex\panel_control.npz` | `35300460da39641568c5e1f31e04e242bd94a94a9166bada3db439efe6764606` |
| `reports_20261027_codex\panel_del_all_feedback.npz` | `da1d8f50bb12278226615fb25d179ce1e306c941b50916b91d2b247f36f96218` |
| `reports_20261027_codex\perturb_w2l.py` | `044bed2b5953149e3719cd8a431c6f4c3867faee3849487636e5908b60d24a36` |
| `reports_20261027_codex\full_basin_results.json` | `8bc2d72bee3c91627bd3534d53c9231cc4276e9782aceeace15fd25318c4327f` |
| `w2l_mujoco\test_w2l_ground.py` | `91df05629c07e48f1a2e62e0d9dd7ac0e82d0cc061a3c66ddcd3ac238cd68d3d` |
| `w2l_mujoco\w2l_ground.xml` | `fa3fd4aa9a6f3fa9b7aaa20c07c6c8a5b83e5d52429b72526cbd7ee4d05b1c59` |

The controller selector is `bilateralrg2`. `recover_w2l.py` checks that saved conductances reach the actual builder, while `finalize_w2l.py` excludes two inert historical keys (`ia_ank_rg`, `ia_ank_pf`). The canonical replay uses the exact-fast neural backend and preserves reference numerical reduction order. These files identify the retained source/configuration chain; the recording package does not preserve a historical compiled MJB or a complete launch-source hash manifest. The current XML hash is therefore not proof of a captured historical compiled binary's identity.

## Scoring definitions and windows

Raw records contain 30,000 samples, **0.000–29.999 s** at 1 ms. No resampling or altered samples were used. The retained metric window begins at `hold + lower + 0.5 = 3 + 1 + 0.5 = 4.5 s`.

- Contact events: heel-force rising crossings above **20 N**, merged if separated by **0.45 s or less**. This is a heel-event counter, not a whole-foot-clearance detector.
- Phase: left primary heel event within consecutive right primary heel events; median across retained event phases.
- Contact phase, event counts, forward COM displacement, contact duty, and mean vertical harness force use **4.5≤t<30 s**.
- Fitted forward COM speed adds the restriction **8≤t<30 s**. The current manuscript's shared 4.5–30 s description must distinguish this extra speed-fit window.
- Combined contact duty is heel **or** toe force > **5 N**.
- Historical supported-stepping gate: finite trajectory, no recorded fall, at least **six primary events on each side**, median phase **0.35–0.65**, and forward displacement > **0.05 m**. It does not itself require the later W2L108 whole-foot clearance/contralateral-loading criteria, or independently test physiological validity.

## Independently recomputed results

All listed scalar metrics agree exactly with their retained per-run JSON values (maximum numerical difference **0.0**).

| Recording | Events L/R | Median phase (cycles) | Forward COM displacement (m), 4.5–30 s | Fitted speed (m/s), 8–30 s |
|---|---:|---:|---:|---:|
| `trial458_d30` | 28/21 | 0.496062992 | 2.246231311 | 0.091785976 |
| `panel_control` | 34/22 | 0.497512438 | 2.238782166 | 0.092370955 |
| `panel_del_all_feedback` | 20/31 | 0.222501060 | 1.995666108 | 0.077014528 |
| `panel_del_Renshaw` | 17/15 | 0.496740548 | 0.080949627 | 0.003089940 |
| `panel_del_commissural` | 25/2 | 0.409516381 | −0.096721799 | −0.003513278 |

Thus the rounded intact/all-feedback phases **0.498/0.223** and displacements **2.239/1.996 m** are correct. All-feedback removal fails the phase gate while retaining cycling and progression. Commissural deletion fails the bilateral event-count and forward-progression requirements; its phase alone is not the failed criterion. Renshaw deletion passes this permissive historical gate despite its small progression, correctly showing that pass/fail alone hides a major functional change.

All **15 numeric trajectory arrays** in each of `repeat_0.npz`, `repeat_1.npz`, `repeat_2.npz`, and `canonical_repeat.npz` equal `trial458_d30.npz` exactly. This independently confirms the four retained replay files, rather than claiming four new simulations. All **12** full saved-scalar basin recordings also pass the recomputed historical gate. `full_basin_w2l.py` applies independent uniform multipliers between 0.99 and 1.01 to numeric saved settings except the fixed 30 s duration; initially zero settings remain zero. These are local parameter perturbations, not physiological validation or global stability evidence.

All listed pathway interventions reproduce their archived pass/fail outcomes: mild/strong extensor-drive reductions, Ia, Ib, Renshaw, crossed-heel, and heel-cascade arms pass; commissural and all-feedback arms fail. Positive `stimN` is subtracted from both extensor RG drive inputs. The retained mild/strong commands therefore reduce those currents by **1.5/3 nA** during **8≤t<9 s**, rather than adding excitatory current. Flexor drive remains unchanged.

## Material protocol qualification: forward pulse is inherited by the control

The retained runner defaults are `pushN=40`, `pushT=8`, and `pushAxis=+x`. Its horizontal force is applied during **8≤t<8.15 s**. `panel_control` omits these keys, and neither `recover_w2l.py` nor `_w2l_wean_eval.py` zeros them. Git's latest modification to this runner is the 27 October campaign commit, not a later manuscript edit. The explicitly named forward-push condition repeats those defaults.

Independent comparison confirms **all 15 numeric arrays in `panel_push_fwd.npz` are exactly identical to `panel_control.npz`**. The backward arm instead reverses the pulse direction. Consequently the forward arm is not a separate perturbation relative to a no-push control. The retained code also gives the nominal control, drive, and deletion arms the common forward pulse. A separate horizontal applied-force trace was not retained, so this protocol conclusion uses the source/default chain and the exact trajectory equality, rather than recovered force telemetry.

Recommended manuscript change: disclose the common source-defined +40 N pulse in the historical fixed-controller comparison and describe the backward arm as reversing its direction. Do not claim a separately demonstrated forward-push perturbation effect against a no-push baseline. The intact/all-feedback coordination contrast remains a matched controller/pathway comparison under the common protocol; the numerical phase/displacement results do not change.

## Support normalization

The legacy scorer divides mean `rigFz` by a fixed **411 N** denominator. `verification.json` records **411.0414 N** for the retained compiled-source XML and explicitly labels the earlier 411 N normalization. Do not silently substitute the latter into the adopted statistic.

For `trial458_d30_ws0.npz`, the independently recomputed mean harness force over **4.5≤t<30 s** is **68.8629698105 N**, giving **16.7549804892%** using the legacy 411 N denominator. The manuscript's historical “approximately 16.76%” is close but is not conventional two-decimal rounding of this particular saved value. Prefer **approximately 16.755%**, with the prior 16.76% preserved in the numerical-change ledger, or retain the existing approximation with the exact recomputed value disclosed in the report. Removing the constant support fraction leaves the vertical PD apparatus active. The separate `trial458_d30` baseline supplies **76.39044084 N**, or **18.58648196%** using the same denominator. These values do not imply harness-free walking.

## Historical IK companion review

`HISTORICAL_IK_REVIEW.md` independently recovers the primary human-model raw recording, original finalist, reference data, scoring implementation, and figure generator. All reported periods, stance fractions, absolute-angle RMSE values, full mean/SD/reference curves, and composite score **−145.0728022062** reproduce exactly. The retained fresh and original recordings share **eight** exactly matching numeric arrays. Two prose corrections are warranted: call the result a **composite kinematic score**, rather than only a cycle-shape score, and state the **2≤t<16 s** analysis window. Historical compiled-model/launch-scorer hash gaps are explicitly recorded there.

## Additional stale status discovered

While locating these sources, `reports_20261027_codex\REPORT.md` was found to open with a later-session correction: the earlier Di Russo native-runtime blocker had been resolved. Its `native_baseline\NATIVE_BASELINE_20261027.md` records a separate author-model replay. The current Appendix's blanket “runtime unavailable” statement is therefore stale and was immediately flagged to the coordinator for a separate primary-evidence check. This does not retrospectively validate the rejected 960-trial MuJoCo adaptation; native source replay and adapted-model locomotion are distinct results.

## Closure after historical correction batch

Read back only the affected current canonical passages; no further source research or generator execution was performed.

- Chapter 6's historical feedback paragraph now states the common source-defined **+40 N forward pulse during 8–8.15 s**, and explicitly says the pair does not establish an unforced baseline. The numerical coordination and progression results are retained under that matched protocol.
- Appendix D's pathway paragraph likewise discloses the inherited pulse and the **exact duplicate forward arm**, while describing the backward arm as a **direction reversal**. It no longer implies independent forward-push-versus-unforced evidence.
- The independent replay paragraph separates contact/displacement analysis **4.5–30 s** from the fitted-speed window **8–30 s**.
- The support paragraph now uses **approximately 16.75%**, explicitly normalized by the legacy **411 N** denominator over **4.5–30 s**. This matches conventional two-decimal rounding of the recomputed **16.7549804892%**; the earlier 16.76% remains preserved above for the numerical ledger.
- The historical human-model comparison now calls **−145.07** a **composite kinematic score**, states the **2–16 s** post-startup window, and says it is not solely a cycle-shape error.

The historical evidence/protocol findings are resolved in prose. Raw recordings and configuration values are unchanged. A minor missing conjunction in the replay paragraph was sent to the coordinator for sentence cleanup; it does not affect any result. Full compilation and rendered-page checks remain separate coordinator verification and are not claimed by this closure.
