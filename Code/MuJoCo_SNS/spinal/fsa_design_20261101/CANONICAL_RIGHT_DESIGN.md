# Explicit right-canonical syn6 candidate

Machine: easteregg2. `mirror_canonical.py` defines a new optional engineering policy. It does not modify the default builder, original K6 files, historical mean-projection code, or published results. No body simulation, Hunt training, or curriculum is performed.

## Meaning of the policy

The candidate defines the left electrical function as the verified-name mirror of the **actual compiled right function**. It does not infer a permutation of independently fitted NMF components or recover a unique biological left synergy basis. The right side retains its existing engineered family classification and phase/time-constant assignments. In particular, it retains the historical mixed-channel interpretation; it does not reclassify the right S5 using a new unilateral stance-fraction calculation.

For every paired neuron, the right capacitance, leak, bias, resting potential, ionic parameters, and numerical integration factors are copied to the left. Synaptic conductance, reversal, and activation-window matrices follow a right-destination ownership rule. A right-target connection is preserved; its left-target counterpart is the exact mirrored connection. If the right connection is zero, a left-only connection is cleared. This is replacement, not averaging.

Shared node parameters and shared-to-shared connections are frozen. For connections into a shared target, the right-side source is canonical and is mirrored to its left-side source. Mixed shared/side connections are included in the edge-change ledger. True unmatched bilateral names and neuron-type or anatomical-muscle-binding mismatches are rejected rather than repaired. Unsuffixed shared nodes are listed explicitly.

Input ports use their actual ordering and connectivity. `BILAT_R/L` spelling is normalized to `BAL_LAT_R/L` for permutation matching; no external port is silently renamed. Ambiguous aliases are rejected. The same target-owned rule is applied to input weights.

Voltages, prior voltages, ionic gates, initialization arrays, time, and phase are not copied or reset. The test initializes unequal side states, applies the policy, and verifies both values and array identities remain unchanged. Parameters can change subsequent dynamics; this is not a claim that an oscillator's later phase is invariant to changed parameters.

Runtime PF family/tau/mixed/DF-channel metadata is explicitly copied from the actual right side. Original `net.syn_W`, source-statistics metadata, and builder graph numerical labels remain historical data. Render or inspect the candidate from its final compiled fields and policy ledger, not the original left NMF labels.

## Correct application order

```python
from fsa_design_20261101.mirror_canonical import (
    apply_canonical_right, default_policy, verify_canonical)
from fsa_design_20261101.graded_profile import (
    apply_profile, default_profile, verify_profile)

net = build_network.build(actuator_names)  # AARL_NET=syn6
apply_canonical_right(net, default_policy())
# Topology changed: construct a NEW accelerated wrapper now.
net.compiled = SparsePairwiseForward(net.compiled)
apply_profile(net, default_profile("syn6"), family="syn6")
verify_canonical(net)
verify_profile(net)
```

`apply_canonical_right(net)` with the default `None` policy is a no-op. A consumable `POLICY.json` is saved with the electrical evidence. The policy requires the dense compiled graph before accelerator construction, because copying/clearing topology can add and remove sparse edges. The later graded profile updates actual accelerated caches and is verified after all transformations.

The inherited `F.tie_weights` operation must make **zero** changes to this candidate. Future half-fit/shared-parameter learning must preserve the declared canonical policy; fitting both sides independently and averaging incompatible PF indices is not authorized by this design. This module contains no optimizer.

## Verified electrical results

- Default syn6 128-step voltage anchor unchanged: `0283aae4d70531aad60a8816611af7c00dc76d7419a646cc8947ccddb12456a4`.
- 794 neurons, with the original population and input-port bindings; no neuron added.
- Canonical projection adds 37 nonzero connections, removes 30 left-only connections, and replaces 122 weights. Every change is named in the ledger.
- Actual neuron/ionic parameters, synaptic conductances/reversals/windows, and input connectivity are exactly bilateral under the verified permutations.
- The inherited averaging operation changes zero pairs, both before and after the graded FSA profile.
- One electrical-step permutation-equivariance error: `2.22e-16`; maximum across 128 steps: `1.78e-15`. A tolerance of `1e-9` was declared. Permuting columns changes floating-point reduction order, so bit-exact permutation equivariance is not claimed.
- The same 128-step candidate voltage trace is bit-exact between dense and accelerated implementations.
- Fresh-process result, source hashes, candidate voltage hash, and all reported checks repeat exactly.
- Original K6 SHA256 remains `63f79c8803b023c9e5ebd5c4c62d571f88d5334df0637058ddc73bfe76c07f6b`; no refit, source mutation, or K8 change.

The final FSA encoder time constant is 3.896229 ms. That remains an explicitly declared per-encoder engineering phase design. Projection equality and correct parameter consumption do not establish body control, a new walking result, full neural stability, or anatomical control coverage. The all-family design gate remains closed.

External evidence: `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\canonical_right_syn6\CANONICAL_RIGHT.json`, `FRESH_REPEAT.json`, and `POLICY.json`.
