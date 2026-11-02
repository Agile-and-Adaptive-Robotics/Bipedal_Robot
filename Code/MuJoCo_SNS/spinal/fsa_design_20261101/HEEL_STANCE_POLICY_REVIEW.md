# Optional heel-routing hypotheses: static design proof

Role: Gait neural pathway composition reviewer. Machine: easteregg2.
No neuronal time steps, body simulation, training, parameter fit, mechanical
change, or default-policy change is included.

`heel_stance_policy.py` provides two separately named, default-off hypotheses:

| Policy | Explicit relay targets per heel | Functional objective |
| --- | --- | --- |
| `A_ipsilateral_extension` | ipsilateral InE positive, ipsilateral InF negative | suppress ipsilateral flexion and disinhibit extension when InF is active |
| `B_contralateral_extension_inhibition` | ipsilateral InE positive, contralateral InF positive | suppress ipsilateral flexion and contralateral extension, matching W2L laterality |

Neither is presented as a uniquely correct interpretation of Syn6. Its original
positive-both-ipsilateral pattern is retained as the historical baseline.
A direct heel-to-contralateral-RG-F excitatory pathway would be a separate
hypothesis; none is added here. Other heel-to-PF pathways remain unchanged.

Public interface:

```python
from fsa_design_20261101.heel_stance_policy import (
    make_policy, apply_policy, verify_policy,
)
apply_policy(net, None)  # exact no-op
policy = make_policy("syn6", "B_contralateral_extension_inhibition")
record = apply_policy(net, policy)  # explicitly opt in, before acceleration
verify_policy(net)  # final scoped compiled-field consumption check
```

Actual builder identity is checked. Missing HEEL/relay cells are rejected, rather
than invented; the default410-neuron S3K therefore cannot silently acquire this
policy. Existing direct heel-to-RG routes are rejected to prevent an undeclared
stacked hypothesis. Applying topology changes after accelerator construction is
rejected; sparse runtime parity has not been claimed in this static audit.

The new operator reuses the actual canonical-right heel-to-InE conductance and
voltage window, then binds the requested targets bilaterally. The reviewed
templates supply0.5 microSiemens,0–5mV source window,+8mV excitation, and−5mV
inhibition. Actual RG inhibitory links have conductance4. No caller gain,
current encoder, capacitance, leak, intrinsic NaP parameter, dynamic state,
geometry, tendon, damping, contact, force law, or physics step changes.
All chemical fields outside the declared heel-to-RG-interneuron scope are
checked byte-exact, and unused profile keys are rejected.

## Static current proof

Proofs use the actual compiled SNS arrays and their clipped conductance formula:

`a=clip((Vpre−elo)/(ehi−elo),0,1)`; `Isyn=g*a*(Erev−Vpost)`.

In the isolated relay fixture, a declared engineering bias holds each IN at
the selected operating voltage without heel current. Other incoming network
paths are excluded. The heel-driven steady IN voltage is calculated exactly
from leak plus synaptic conductance, followed by the current change through its
actual inhibitory RG projection. This is a local steady operator, not a claim
about the equilibrium or rhythm of the complete recurrent network.

At a2.5mV IN/RG operating point and0.75nA right heel input, all three reviewed
families produce the following signed changes:

| Hypothesis | Right RG-F | Right RG-E | Left RG-F | Left RG-E |
| --- | ---: | ---: | ---: | ---: |
| A | −2.302325581nA | +3.139534884nA | 0 | 0 |
| B | −2.302325581nA | 0 | 0 | −2.302325581nA |

The left-heel operator mirrors these targets. At zero encoded heel current the
isolated effect is zero. Grid checks additionally cover0,0.125,0.75,1.5,5,and20nA
at0,2.5,and4mV operating points. Signs and laterality do not reverse;5 and20nA
give identical saturated pathway currents. A cannot disinhibit an already-silent
InF projection: at the0mV relay operating point its extension-current change is
zero. This condition is recorded rather than treated as a failed sign.

The0–0.75nA interval comes from the actual saved W2L heel caller. Larger currents
are engineering saturation probes, not admissible inputs for that saved caller.
Syn6's saved standing heel gain remains zero; the0.75nA check is a separate
expression fixture and does not silently enable its physical input.

For S3K the expression fixture explicitly enables source `full_rules=1` and
`heel_rge=0.5`. It inherits no standing result. The default no-HEEL constructor
was separately checked to reject an enabled policy, while `None` remains exact.

## Pulse versus held contact

The policy changes graph targets only. It does not detect heel-strike onset.
An explicit upstream event encoder is required for pulse-only reset use.
A held high heel input continues to apply the high static operator throughout
contact; it is tonic bias, not a one-time reset. Low-input steady current is
zero, but membrane transients and phase-transition timing are not tested here.
No complete-network, sparse-timestep, locomotion, or stability result is inferred.

Evidence: three compact JSON files in
`D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\gait_full_composition_v1\heel_AB_static`.
Each includes before/after parameter hashes, exact changed fields, final scoped
verification, historical-baseline operators, both hypotheses, domain grid,
source hashes, and explicit zero neural/physics step counts. Central Ia/II/Ib
rule additions remain a separate, unapplied design-review matrix in
`GAIT_TOPOLOGY_SOURCE_REVIEW.md`.
