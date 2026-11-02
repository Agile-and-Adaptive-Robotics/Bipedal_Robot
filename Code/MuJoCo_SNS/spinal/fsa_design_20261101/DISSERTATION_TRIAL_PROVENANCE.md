# Were central feedback gains zero during the dissertation trials?

Session 1. Read-only trace; no simulation rerun, training, or manuscript edit.

**Not for the full-s3k baseline/ablation trials.** Earlier conclusions about the later standing export cannot be generalized to the submitted results.

## Specific result traced

Archived submitted Appendix D reports full-s3k −186.4 and syn6 −200.4. `Code\MuJoCo_SNS\spinal\prune_analysis.md` records the corresponding unrounded WALK baselines −186.41420708638907 and −200.37955834081458. The recorded raw results are `prune_results_s3k.jsonl` and `prune_results_syn6.jsonl`, selecting `config=full`, `mode=WALK`.

`prune_matrix.py` explicitly loads each family's own curriculum configuration before applying any ablation. For full s3k it merges v10 base multipliers with `reports_20260923/s3k_trial34_full_params.json`, sets Renshaw to 0.5, and calls `_curriculum.set_stage(5, WALK_P)`. That setter assigns the central gains before the runner constructs its network.

| Full-s3k saved parameter | Value |
| --- | ---: |
| ia_f_central | 0.3985916558720821 |
| ii_f_central | 0.38389852925373125 |
| ii_e_central | 0.09300453229764793 |
| ib_e_central | 0.24354654619098992 |
| ia_in | 0.8250460785281135 |
| full_rules | 1.0 |
| renshaw, explicitly set by campaign/setter | 0.5 |
| heel_rge | 0.37530421840765804 |
| toe_rge | 0.14572281192093522 |

The separate Ia/II/Ib ablations then zero selected gains intentionally. The full condition must not be conflated with an ablation or the later 410-neuron standing condition.

## Historical syn6

The same campaign loads `curriculum_syn6_stage5.json` through `_curriculum_syn6.set_stage(5, WALK_P)`. Its saved heel and toe gains are 0.12312339160142002 and 0.10178084469071401; ib_rge is zero. This loader enables contact drive and explicitly sets the pooled Ib central pathway from that zero value. It does not enable the stock central Ia/II projections. Thus the heel graph is not merely a silent-standing issue for this historical recipe. Renshaw is explicitly disabled in the syn6 setter.

These are source/configuration/result provenance checks. They are not a newly rerun historical neural trace or a reconstruction using a guaranteed byte-identical old source revision. Historical-code reconstruction must be labeled separately. Other dissertation figures and distinct later comparisons require their own recipe trace rather than inheriting this conclusion.

## Why the drawings differed

The supplied EB475WS4 JSONs point to later October 31 standing checkpoints, not the above pruning campaign. Their saved full_rules and central gains are zero. The s3k 410-neuron drawing therefore cannot illustrate the historical full-s3k dissertation trial. Files and historical results are preserved; diagrams must identify which actual condition they represent.
