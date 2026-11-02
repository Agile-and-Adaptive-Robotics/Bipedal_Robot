# Existing central afferent current budgets

Machine: **easteregg2**. Helper: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\central_afferent_budget.py`.

This is a new optional electrical contract, not a missing-count correction or a trained walker. The existing w2lvar builder already divides its central Ia/II gains by `_nF` and Ib gain by `_nE`. Its 216 direct Ib central edges do not imply 216 independent full-strength copies.

The helper reconstructs each saved biological standing circuit through the actual runner's nested builder, before any body binding. The three builders now have authorized default-off profile hooks, so their file hashes differ from the saved runs. Before proceeding, the helper verifies all five archived compiled matrix hashes, saved bilateral projection, neuron parameters, actual membrane/ionic fields, timestep, and input mappings against the frozen graph exports. Any numerical mismatch stops the check. This lineage change is recorded explicitly.

## Input envelope and function

Canonical assumptions are gamma coactivation0, vestibular proprioceptive gain0, and muscle-afferent scaling1. Actual runner/AFF/MOD values give conservative external-current bounds:

- Ia: `1.6*(1+0.8)=2.88 nA`.
- II: `1+1.2*(1+0.6)=2.92 nA`.
- Ib: `0.5*1.5*1.5=1.125 nA`.

These follow clipped velocity/length in[-1,1], force/Fmax in[0,1.5], and drive/stance in[0,1]. They are uniform envelopes, not simultaneous physically attainable maxima. In particular, some nonextensor II channels have a smaller bound because their variable gain is halved. Changes to gamma, afferent scaling, or other input-generation rules require a new declared domain; nothing silently clamps production inputs.

Only existing direct Ia/II/Ib→RG/PF/lamination edges are selected. An isolated passive encoder's maximum voltage follows its actual input coefficient and leak; the actual presynaptic threshold/saturation then determines its maximum conductance fraction. Existing within-group relative weights are preserved. Each group is one postsynaptic cell, one afferent type, and one current sign at Vm2.5mV. The explicit engineering starting budget is0.5nA in magnitude per group. This is an incremental NaP input-current contract, not a steady-voltage gain target. Budgets across types add; Ia0.5 plus II0.5 can give1nA to one target. Actual current changes with the target's membrane voltage.

One verified right half determines its mirror. Names, full source/target membrane and ionic parameters, reversals, windows, source bounds, weights, and complete group membership must match. Unmatched groups remain unchanged and are listed. HUNT_AF paths are excluded. Composite INs, including LBIN, are inventoried separately and are not assigned a direct-afferent bound automatically. None of these saved standing circuits has an active LBIN→central path requiring such qualification.

## Results and scope

| Saved circuit | Existing selected edges | Post/type/sign groups | Scope |
| --- | ---: | ---: | --- |
| normalized w2lvar |520|24|New current-budget electrical checks pass; no body result.|
| s3k |0|0|Inventory/default-preservation only; no absent pathway invented.|
| syn6 |0|0|Inventory/default-preservation only; no absent pathway invented.|

In w2lvar, the uniform-envelope old reference bounds at RG Vm2.5mV are Ia1.67904nA and II1.70236nA into RG-F, and Ib0.61875nA into RG-E. These are algebraic envelope bounds, not recorded standing currents. The new per-type0.5nA targets multiply existing weights by approximately0.297789,0.293710, and0.808081, respectively. The new function intentionally differs from the previous aggregate conductance contract.

Native sparse on/off comparisons at0/25/50/100% of the declared source envelope reproduce the summed-current calculations within1.24e-14nA. Across120 actual base-network steps, selected passive encoder voltages stay inside their derived bounds, the neural states remain finite, and cached/dense conductances remain consistent. All nonselected weights, membrane fields, ionic parameters, and timing constants are preserved. The final w2lvar result/profile/ledger repeat byte-exact in a fresh process.

The s3k/syn6 `electrical_passed` bookkeeping includes successful inventory and exact default-preservation checks; it is **not evidence of a tested central function where zero edges exist**. Every record leaves global family verification false.

## Application and evidence

For an already constructed compatible base network:

```python
import central_afferent_budget as CB
CB.install(net, CB.default_profile())
CB.verify_runtime(net)
```

`CB.install(net,None)` is an exact no-op. This is a compiled overlay; rebuilding the network requires explicit reapplication. It updates both actual dense and sparse fields. Networks with Hunt's dynamic `apply_weights` policy are rejected: no Hunt runtime integration is claimed. Production callers must preserve the stated physical input-domain assumptions and check any additional masks.

Final external records:

- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\central_afferent_budget\w2lvar\domain_verified\RESULTS.json`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\central_afferent_budget\w2lvar\domain_verified_fresh\RESULTS.json`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\central_afferent_budget\s3k\domain_verified\RESULTS.json`
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\central_afferent_budget\syn6\domain_verified\RESULTS.json`

Each folder also contains the explicit PROFILE.json and applied ledger. Earlier compact checks are retained as lineage. No source models, accepted winners, manuscript files, or ZCode experiments changed. Physics, body binding, Hunt fitting, and curriculum execution are all zero; training gates remain closed.
