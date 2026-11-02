# Dissertation Gait recipes: topology and static FSA review

Role: Gait neural pathway composition reviewer. Machine: easteregg2.

The dissertation prune trials must be distinguished from later standing exports.
In particular, the dissertation S3K recipe enables central sensory feedback and
explicit reflex interneurons. The 410-neuron standing/default audit does not
describe that recipe. This review gives separate full, combined-prune, and
pruned-retune conditions. No historical result or parameter file was changed.

## Source recipes and observed historical metrics

| Condition | Current reconstructed neurons/nonzero pairs | Direct central Ia / II / Ib pairs | Historical record |
| --- | --- | --- | --- |
| S3K full | 878 / 8340 | 152 / 368 / 216 | full/WALK kine score−186.41420708638907 |
| W2L full | 888 / 7868 | 152 / 152 / 216 | full/WALK kine score−144.46440311078322 |
| Syn6 full | 794 / 3037 | 0 / 0 / 0 | full/WALK kine score−200.37955834081458 |
| S3K combined prune | 868 / 8084 | 152 / 368 / 0 | combo/WALK kine score−176.08919108631918 |
| S3K pruned walk1 fix | 868 / 8084 | 152 / 368 / 0 | saved optimization score−187.00327631586913 |
| S3K pruned walk2 fix | 868 / 8084 | 152 / 368 / 0 | saved optimization score−197.51274884782532 |

The first four metrics come from original `prune_results_{family}.jsonl`
records. The last two come from `curriculum_pruned_walk_{1,2}_fix_20261003.json`;
their optimization objectives are not silently equated with a raw WALK metric.
These values were read, not recomputed. Frozen/one-sided behavior in the saved
metrics is retained; a filename containing WALK does not establish walking.

`gait_dissertation_topology_audit.py` executes only source-extracted
`optuna_walk.load_fitted_baseline`, `set_params`, and the relevant curriculum
`set_stage(5,p)` function. It merges the actual v10 multipliers with the actual
family recipe, following `prune_matrix.py:67–134`. The S3K source recipe is
`reports_20260923/s3k_trial34_full_params.json`, not a later standing winner.
The combined prune follows literal `prune_combo.COMBO_KEYS` and
`--no-interleg`. Pruned retunes execute the actual extracted
`_curriculum_pruned._config`, validating saved PIN/prune fields against source.

The original JSONL records do not supply historical source/model/config hashes.
Consequently these are **current-code static reconstructions of the saved
recipes**, not verified historical byte-identical networks or reproduced body
results. Current source, model, recipe, and result hashes are saved. The current
runner default model is parsed as XML for actuator names; no MuJoCo model is
loaded. `connectome_gains.json` is absent now; historical external overlays and
environment overrides are not proven. Runner/main is neither imported nor called.

## Interneurons and source flags

Full S3K and all three pruned S3K conditions contain 92 each of IaIN, RC, IIX,
IIIN, and IBIN, plus 10 IBEXC cells. Full S3K has `full_rules=1`,
`ia_in=0.8250460785281135`, and `renshaw=0.5`. `ia_in>0` enables its topology;
that value is not itself every IaIN synaptic conductance. The pruned walk2
retune has `ia_in=0.3711788870663949`, still enabling the same population set.

W2L contains the same five 92-cell relay sets and 10 IBEXC cells through its own
builder. Syn6 has 92 each IaIN, IIX, IIIN, and IBIN, plus 10 IBEXC, but no RC.
Its `ia_in=0` and `full_rules=0` flags do not remove its independently implemented
relay populations. Names of shared flags are insufficient to infer topology.

S3K full central conductances are Ia-F0.3985916558720821,
II-F0.38389852925373125, II-E0.09300453229764793, and
Ib-E0.24354654619098992 microSiemens per actual source edge. These paths are
active in the reconstructed dissertation recipe. Combined pruning removes the
direct central Ib path and contact/commissural features, while retaining the
central Ia/II routes, local autogenic Ib circuitry, and IBEXC pools. It is not
global removal of all Ib feedback.

Syn6's saved stage5 caller has heel gain0.12312339160142002 and toe gain
0.10178084469071401; `ib_rge=0`. Its contacts are therefore not equivalent to
the later zero-contact standing caller. Its direct central Ia/II/Ib projections
remain absent in this actual recipe, while local afferent/reflex routes remain.

## Actual heel connections

Full S3K has selective ipsilateral heel routing: heel excites InE and inhibits
InF. Its final heel-to-InE conductance is0.24354654619098992, whereas
heel-to-InF inhibitory conductance is0.37530421840765804. The difference is
observable in final matrices: `_build_pf` adds a later heel-to-InE definition
using `ib_e_central` (`build_network.py:602–609`), overwriting the earlier
`heel_rge` definition in pinned SNS's one-entry-per-pair compiler. Do not report
both final connections as having the configured heel_rge conductance.

W2L heel excites ipsilateral InE and contralateral InF. Syn6 heel excites both
ipsilateral InE and InF. These are different functional hypotheses, not a rule
that both models share. Existing optional A/B heel policies were **not applied**
to these historical reconstructions or to the FSA candidates below.

## Static opt-in FSA compatibility

Each of the six reconstructed networks was cloned into a separately identified
NEW optional graded FSA candidate. Each passes:

- 368 signed terminal-pool equilibrium checks; largest error
  4.440892098500626e−16mV against the declared engineering contracts.
- 276 actual muscle-afferent port bounds derived from current source formulas:
  Ia0–2.88, II0–2.92, and Ib0–1.125nA, with gamma0 and vest_prop0.
- Actual input-adapter upper-endpoint identity, pA-to-nA conversion, and rejection
  of a current beyond each afferent type's declared bound.
- Final sparse parameter-cache identity with the applied dense profile.
- Unchanged NaP parameters and nonencoder RG/PF/IN/MN capacitances, plus an
  unchanged historical-baseline parameter hash.

The source-derived encoder candidate has tau 3.8962288367 ms. At a 2 ms update,
the actual discrete transfer function has 10 degrees phase lag at the declared
7.155635 Hz engineering band. This is an encoder design allocation, not a new
whole-network clock or proof that a whole closed loop meets that phase budget.
NaP bursting, PF shape, and source intrinsic periods are preserved.

The terminal voltage contracts are isolated engineering targets, not newly
observed physiological responses. The bounded signal checks cover muscle
ports, not the entire descending, contact, event, posture, or balance caller.
There were zero neuronal time steps, zero body steps, zero training trials,
and no dynamic sparse/full-network parity or stability claim in this audit.

## Central-current inventory is separate and unapplied

The existing current-budget helper was used to propose/inventory matched
existing routes, without installing new budgets or edges. Full S3K has 32
qualifying direct groups/736 edges; combined-pruned S3K retains 24 groups/520
edges; W2L has 24 groups/520 edges; Syn6 has no qualifying direct groups.
Composite LBIN/IN paths are listed separately and are not silently treated as
direct afferents or certified by these terminal checks.

For example, full S3K's 19 Ia sources into right RG-F yield a 23.9920289503 nA
upper-reference current at the helper's declared 2.5 mV operating point; its
19 II sources yield 23.4285594433 nA. These are mathematical source-envelope
evaluations, not measured simultaneous physical currents or a demonstrated
cause of failure. The proposed 0.5 nA per-target/type/sign budget remains an
engineering proposal, not an applied physiological correction.

No canonical mirror, central-current budget, Ib coincidence, heel policy,
tendon compliance, damping repair, contact change, Hunt fit, or body curriculum
was applied in this checkpoint. Biological route completeness and new functional
rules remain separate from the passing unit/scaling checks. Training release
remains closed pending Ben's explicit approval.

Final evidence is the six `{family}_{condition}_DISSERTATION_PRUNE_TOPOLOGY.json`
files in `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\gait_dissertation_prune_topology`.
The older three files without a condition suffix are preserved initial
topology-only checkpoints; use the six condition-specific files for this review.
