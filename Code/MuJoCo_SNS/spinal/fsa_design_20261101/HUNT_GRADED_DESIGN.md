# Hunt and four-cell local-reflex design

Session 1, easteregg2. This is an opt-in electrical design, separate from the saved standing/stepping conditions and training records.

## Runtime finding

The saved Hunt/PT stepping condition has 1033 neurons; the earlier four-cell dynamic condition has 1027. In `hunt_training/controller.py`, `install()` clears all MN input conductances and classifies a selected set of fitted reflex edges. `apply_weights()` restores that parameter set on every command update. Homonymous IIX-to-MN and antagonistic IIIN-to-MN edges are allocated in the builder, but omitted from that restored set.

For the saved Hunt/PT condition, all 92 homonymous IIX-to-MN connections therefore have zero conductance at command values 0, 0.5, and 1. This does not mean that the separate HUNT_AF timing features contain no length information. The finding concerns the local II-to-MN pathway, not all possible type-II-related input to the network.

Static audit: `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\HUNT_PT_RUNTIME_PATHWAYS.json`. It excludes the separately body-calibrated PT output columns, whose gain writes are disjoint from these native afferent classes.

## Opt-in correction

`hunt_graded_profile.install(net, profile)` starts the four declared local sensory functions from the existing canonical template, derives the signed-pool conductances, and updates the **authoritative runtime** `reflex_gain` and `reflex_edges` arrays. This restores the 92 existing homonymous excitatory II paths and prevents the next command update from clearing them. It also initializes the Ia homonymous, Ia reciprocal, and Ib inhibitory pools to the declared FSA functions. Historical fitted weights remain untouched on disk.

It introduces no neuron or new anatomical connection, preserves the NaP/PF/HUNT_AF and geometry policies, and leaves the conflicting antagonistic-II policy inactive. The latter is a separate biological/design decision, not silently included in this correction. A `None` profile is an exact no-op.

The candidate uses the graded profile's stated electrical engineering contracts: a 1 nA calibration pulse, a signed 0.5 mV response about a 2.5 mV MN operating point, and the discrete-time fast-encoder setting. These are testable initialization targets, not physiological measurements. Local gain changes are a new controller condition to be trained later, not a retrofit to old results.

## Actual electrical checks

Both saved source conditions were constructed in fresh processes, with no body integration:

| Condition | Neurons | Old active homonymous II paths | New active paths | Soleus MN response, R/L |
| --- | ---: | ---: | ---: | --- |
| Hunt/PT saved stepping condition | 1033 | 0 | 92 | +0.50000000000008 mV / same |
| Four-cell dynamic RF/HAM condition | 1027 | 0 | 92 | +0.50000000000008 mV / same |

The local reflex parameters remain present after command updates 0, 0.5, 1, and back to 0. All 356 modified left/right edge pairs are checked for equal absolute conductance. In the actual compiled networks, the isolated soleus II-to-IIX-to-MN trajectories match between sides exactly. The 2.5 mV baseline and 0.5 mV response agree with the declared endpoint, and the default-None configuration remains exact.

These checks isolate the two relevant native pathways in the real compiled circuit. They do not replace it with a two-neuron toy model. PT mechanical calibration, intrinsic rhythm, other central/biarticular interactions, and the full control-coverage gate remain separate.

Evidence:

- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_local_reflex\Hunt_PT\RESULTS.json`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_local_reflex\Hunt_PT\APPLIED_PROFILE.json`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_local_reflex\fourcell_RF_HAM\RESULTS.json`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_local_reflex\fourcell_RF_HAM\APPLIED_PROFILE.json`

## Use and later training

Construct the actual saved Hunt network and finish its source weight installation first. Then call the opt-in installer with the explicit graded profile. The correction acts through the existing runtime gain buffers, rather than an after-the-fact MN current clamp. Call `verify_runtime(net)` after an intentional command/geometry update.

Later Hunt refinement must parameterize the declared functional gains and save its resulting profile as well as the fitted circuit weights. Do not train arbitrary raw reflex weights and then reload this initializer over them: that would discard the learned local functions. New training/replay integration must explicitly distinguish initialization from consumption of a trained profile. No such training has started; the all-family gate remains closed.
