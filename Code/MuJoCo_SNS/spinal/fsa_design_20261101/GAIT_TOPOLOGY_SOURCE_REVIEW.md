# Gait neural pathway composition: source topology audit

Role: Gait neural pathway composition reviewer. Machine: easteregg2.
This audit constructs neural graphs but advances neither neurons nor physics.
No training, source-wiring edit, mechanical change, or release approval occurred.

Scope: the raw-default and later standing graphs below are not the dissertation
prune-trial configurations. In particular, dissertation S3K enables central
feedback and explicit interneurons. Use `GAIT_DISSERTATION_RECIPE_REVIEW.md` for
the separately reconstructed full, combined-prune, and pruned-retune recipes.

## Actual interneuron lamination

Counts below are compiled nonzero pairs, rather than connection definitions that
may contain zero conductances or duplicate source/target entries.

| Constructor | Neurons | Nonzero pairs | IaIN | RC | IIX | IIIN | IBIN | IBEXC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| S3K raw default | 410 | 1186 | 0 | 0 | 0 | 0 | 0 | 10 |
| S3K requested `full_rules=1` only | 686 | 2334 | 0 | 0 | 92 | 92 | 92 | 10 |
| S3K archived standing | 410 | 836 | 0 | 0 | 0 | 0 | 0 | 10 |
| W2L raw/default | 888 | 7864 | 92 | 92 | 92 | 92 | 92 | 10 |
| W2L archived standing | 888 | 7570 | 92 | 92 | 92 | 92 | 92 | 10 |
| Syn6 raw/default | 794 | 3033 | 92 | 0 | 92 | 92 | 92 | 10 |
| Syn6 archived standing | 794 | 2966 | 92 | 0 | 92 | 92 | 92 | 10 |

The 410-neuron S3K graph is a reduced direct-path implementation, not a complete
interneuron-laminated reflex model. It has 92 excitatory Ia-to-MN pairs, 436
inhibitory Ia-to-antagonist-MN pairs, 92 direct excitatory II-to-MN pairs, and 92
direct inhibitory Ib-to-MN pairs. W2L and Syn6 retain the 92 monosynaptic
excitatory Ia-to-MN connections; their II and autogenic Ib paths use relays.

`build_network.py:257,272` independently gates RC and IaIN creation using
`renshaw>0` and `ia_in>0`. Both default to zero (`params.py:119,218`).
`full_rules=1` alone changes II/IIIN/IBIN motifs (`build_network.py:760–791,
890–920`), without enabling IaIN or RC. The independent IaIN path and its
PF/RC modulation are at lines 872–915. Any fully laminated S3K candidate must
explicitly declare these additional enables, its conductances, and the resulting
new graph; no such candidate was silently selected here.

W2L's own constructor installs its relay motifs and resets `full_rules` to zero
through `VARIANT_G` (`build_network_w2lvar.py:128–143,701+`). Asking this builder
for `full_rules=1` therefore has effective value zero and does not change its
888-neuron graph. Syn6's own motifs are likewise independent of this flag;
setting it to one leaves 794 neurons and the actual connections unchanged. Its
current graph has no RC populations. A shared flag name does not establish
identical implementation across families.

## Heel pathways: laterality matters as well as sign

All following first-hop conductances are 0.5 microSiemens with excitatory
reversal +8 mV unless stated otherwise. The RG relay links are inhibitory,
conductance 4, reversal −5 mV (`E_REV_INH=-E_HI`): InE inhibits RG-F;
InF inhibits RG-E.

| Source policy | Heel to relay | Two-hop effect in the active synaptic window |
| --- | --- | --- |
| S3K `full_rules=1`, with contacts separately enabled | ipsilateral InE positive; ipsilateral InF negative | inhibits ipsilateral RG-F and disinhibits ipsilateral RG-E |
| W2L actual source/default | ipsilateral InE positive; **contralateral** InF positive | inhibits ipsilateral RG-F and contralateral RG-E |
| Syn6 actual source/default | ipsilateral InE positive **and ipsilateral InF positive** | inhibits both ipsilateral RG half-centers |

Source locations: S3K `build_network.py:512–534`, W2L
`build_network_w2lvar.py:279–286,338–344,378`, Syn6
`build_network_syn6.py:433–438,568–574`. W2L's comment describing inhibition of
the partner flexor conflicts with its actual InF-to-RG-E target. The source
connections, not that comment, determine the functional sign.

This does not establish a unique Syn6 correction. Separate optional hypotheses
are A: selective ipsilateral extension, matching the S3K full-rules pattern;
and B: ipsilateral stance/contralateral-extension inhibition, matching W2L.
A direct heel-to-contralateral-RG-F excitatory pathway is a third hypothesis,
not silently stacked with either. Changes require isolated consumed-matrix
proofs before any later user-authorized body evaluation.

Syn6's saved caller has heel/toe gain zero, despite the edges existing. W2L's
saved heel/toe currents range 0–0.75 nA. A silent source port does not validate
the biological function of its stored connections. Held heel input continues
to bias these pathways; none of these connections alone implements a heel-strike
onset pulse. Pulse and held-contact functions must be declared separately.

S3K can add an excitatory heel-to-InF definition after an inhibitory definition
(`build_network.py:532–534`). Pinned SNS `networks.py:465` appends connection
definitions, while `compilers.py:330–333` assigns one matrix element per pair.
Such duplicate definitions can overwrite rather than sum opposite-reversal
channels. Audits must inspect final matrices, not count source declarations.

## Muscle afferents: local versus central, default versus optional

| Family/default and archived standing | Direct central Ia | Direct central II | Direct central Ib |
| --- | ---: | ---: | ---: |
| S3K | 0 | 0 | 0 |
| W2L | 152 | 152 | 216 |
| Syn6 | 0 | 0 | 0 |

W2L's central routes are present in both raw and selected graphs, independently
of the standing PF-to-MN output mask. Flexor-family Ia/II target the same-side
RG-F, InF, hip-F PF, and knee-F PF. Each edge has conductance 0.53/19 =
0.027894736842105264. Extensor-family Ib targets RG-E, InE, hip-E PF, and knee-E
PF; each edge has 0.5/27 = 0.018518518518518517. This is count-normalized source
wiring (`build_network_w2lvar.py:595–608`), not an absence of central feedback.
It does not implement extensor II-to-E in that block or joint-specific ankle,
lumbar, frontal, or transverse targeting.

S3K contains optional F Ia/II-to-F and E Ib/II-to-E central routing at
`build_network.py:805–850`. Corresponding gains all default to zero
(`params.py:275–278`); they are independent of `full_rules`.
Targets are RG, RG interneuron, and own joint PF when configured, with phase-PF
fallback otherwise. Optional relayed PF routes appear at lines 927 onward.

Syn6 implements optional extensor Ib-to-E central routing
(`build_network_syn6.py:765–778`) and optional Ib-to-LBIN-to-E routes
(lines 590–599,737–740), both default off. Its own muscle wiring does **not**
implement the S3K F Ia/II-to-F or E II-to-E central blocks. Setting a similarly
named unused gain or `full_rules` is not proof of a consumed implementation.
Its local afferent/MN/IN routes remain present.

Deafferentation is a caller condition, not necessarily a graph deletion.
Runner `use_aff=False` skips raw muscle afferent encoding; separate contact
channels require separate inspection. `AARL_DEAFF` also masks selected contact
inputs. Hunt removes inherited sensory-timing projections and supplies its
separate HUNT_AF channels; it retains local sensory paths. Do not describe all
of these conditions as the same full/deafferented network.

## Proposed typed central rule matrix: review only, not applied

| Source class | Candidate central targets | Existing precedent | Binding needed before any new wiring |
| --- | --- | --- | --- |
| Functional flexor Ia, II | own-side F RG/IN and relevant F PF | S3K optional; W2L active | source-length/velocity/current adapter and signed contribution |
| Functional extensor Ib | own-side E RG/IN and relevant E PF | S3K/Syn6 optional; W2L active | actual load encoding, phase gate, and summed-current budget |
| Functional extensor II | own-side E RG/IN and relevant E PF | S3K optional | source joint/function targets and gain/current contract |
| Biarticular muscle | targets matching each signed joint action and phase | four-cell RF/HAM primitive design | do not classify RF/hamstrings solely by one muscle-family name |
| Abductor/adductor, rotator, lumbar, foot | appropriate task/axis PF plus explicit sensory function | coverage inventory only | source moment-arm signs across poses, role, and closed-loop objective |

No blanket all-to-all projection or arbitrary new gain is justified by this
table. Preserve each muscle's source contribution distribution, then budget
summed conductance/current per functional destination. K6 channels can mix
joint actions; primary muscle names or same-index bilateral channels do not
prove functional equivalence. The canonical-right candidate is an explicitly
engineered mirror policy, not recovery of a unique bilateral NMF alignment.

## Evidence and interpretation

`gait_topology_audit.py` produced three compact JSON files under
`D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\gait_full_composition_v1\topology_source_audit`.
All three selected graph hashes exactly match the EB475WS4 files supplied for
review. JSON includes full source hashes, actual masks, raw/default and requested
full-rules matrices, named heel paths, direct sensory counts, and source globals.
Neural steps = 0; physics steps = 0.

Companion caller-domain supplements in `caller_domain_supplements` cover all
376 S3K or 382 W2L/Syn6 ports using actual archived input matrices, static
activation targets, saved scales, and runner formulas. Source/endpoint checks
do not imply that every pathway has been activated in a body simulation.

Earlier unit/FSA/scaling proofs establish operator consumption, units, saturation,
and bounded electrical behavior. They do not prove that the original graph has
all intended biological interneurons, correct heel laterality, or complete
central sensory rules. Those claims are explicitly separated by this audit.
