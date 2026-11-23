# Human Gait2392 stance and joint-angle revision

## The original score

The 280-trial feedback search maximized the sum implemented in spinal/kine_ref.py and spinal/wean_tune_20261016.py. Its waveform-error term first centered the model and recorded joint curves separately. It also rewarded joint excursion and intra-cycle phase agreement, and penalized mean-posture error, knee minimum, period error and variability, stance duration, simultaneous loading, and between-leg phase. The higher bilaterality score required at least three cycles per leg. The objective did not minimize the mean of the six unshifted joint RMSE values. Reconstructing each component from the archived data exactly reproduced both composite scores. Matching excursion reduced the penalty by 50.41 points, whereas matching stance reduced it by 0.57. The composite improved by 44.66 points. Under the original convention the right/left ankle RMSE increased from 6.82/10.62° to 19.40/27.19°.

## Physical force and affected results

In the human Gait2392 runner each foot contact delivered the same normal force to the heel and toe sensor channels. The sum saved as “contact” counted it twice. A fresh replay measured all four calcaneus/toe–ground contact pairs separately. Across all 8,000 sampled steps the maximum difference between the independently measured upward force and half the recorded signal was 0 N. The feedback entering the old neural controller and the archived trajectories are preserved. Evaluating the human recording and simulation at a physical 50 N gives the following results.

| Human controller in a 16 s trial | Mean six-joint RMSE | Right/left stance | Right/left cycles |
| --- | ---: | ---: | ---: |
| Before feedback retuning | 17.78° | 76.2/66.9% | 11/10 |
| First weighted-score search | 21.54° | 50.5/44.8% | 12/13 |

The older 17.78° and 23.36° values used physical 25 N for the simulated foot and 50 N for the recorded human foot. Every numerical change above is a reanalysis of the same motion. An independent repeat also verified the matched raw trajectory hash. Changing the cycle boundaries changes the averaged joint-angle curves, even though the movement is the same.

Appendix D now reports corrected force thresholds, gait periods, joint-angle table and overlay, loading-cycle counts, stance fractions in the short and long trials, and robustness to the original eight parameter perturbations. At the defined 40–75% criterion, three of the original eleven supplementary-drive settings qualify. Five of eight perturbations of the weighted-score winner qualify. The 22 s weighted-score trial falls outside this range on the right. The separate spiking/nonspiking comparison and the external-reference results were also recomputed at their stated physical thresholds. The engineering record **STANCE_AND_TABLES_EVIDENCE.json** gives both sets of numbers, including per-joint, drive-scale, perturbation, and reference values. The CAD, robotbody, and W2L studies use separate recordings and scoring paths.

## Typesetting

Tables 4.2, 4.3, C8, and C9 produced all eleven underfull horizontal-box messages in the previous full local LaTeX log. Only the alignment of the wrapping text columns changed in these four tables. Their captions, symbols, entries, and numbers agree with the source. The table-only local and live builds both cleared the diagnostics. The final local print check covers the whole dissertation, including the updated gait figure. The output-file checks are in the evidence record.

## Editorial check

The advisor's question “Why is this important? What did you do about it?” applies directly to the former account of two RMSE values and stance percentages. Chapter 6 now identifies the trade between joint excursion and ankle posture. Appendix D identifies the objective term responsible, defines how stance is measured, and distinguishes the two thresholds. Chapter 7 identifies an independent test on other recordings.

## Three-objective controller fit

The corrected 50 N definition governed 531 new evaluations. Of those, 428 passed the stability and bilateral cycle criteria. Zero met both the 17.78° and 58–64% targets. The lowest angle error compatible with bilateral stance was 18.48°, with 62.3/61.4% right/left stance in an independent 16 s replay. In constant drive, those fractions were 77.5/75.4%. Nominal 22 s and independent 28 s tests of the same settings gave 20.89° and 63.2/53.7%, and 21.20° and 69.3/62.3%, respectively. The candidate with the lowest overall angle error (16.98°) gave 80.8/62.7% right/left stance at 16 s.

The longer trials and parameter perturbations were withheld from selection. Appendix D compares the six joint errors and plots the three distinct trajectories. The archive in Code/MuJoCo_SNS/spinal/physics_fidelity/harness_sensitivity/pareto_reference contains individual trial settings, measured loading, results, the Pareto front, and independent validation.
