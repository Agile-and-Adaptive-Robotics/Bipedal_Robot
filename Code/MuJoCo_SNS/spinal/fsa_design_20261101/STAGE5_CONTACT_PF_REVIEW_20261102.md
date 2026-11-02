# Stage5 actual contact-PF route review

Session 1, easteregg2. Read-only actual construction from the **earlier `pt_low_s4.json` configuration**, not the later ground9 final output profile. No neural or body stepping, parameter fitting, or implementation change was performed by this audit.

The separate training record defines stage5 as **ground contact-PF refinement with a measured support-reduction ladder**. It is not a renamed perturbation stage. Actual `physical_eval.build` restores existing HEEL/TOEDF→PF slots from `original_g*contact_pf_scale`, before FSA composition. Exactly six pairs change; all other conductances match stage4 and the six persist through commands0/0.5/1/0. [Exact source/configuration/compiled proof](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/STAGE5_CONTACT_PF_AUDIT.json>), [audit code](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/audit_stage5_contact_pf.py).

| Actual per-side restored route | Functional effect | Gain at scale 0.25 /1, µS |
| --- | --- | --- |
| HEEL→PF_IN_HIP-E | Excites the extensor-labelled inhibitory relay; relay inhibits PF_HIP-F | 0.125 /0.5 |
| HEEL→PF_IN_KNEE-E | Excites the extensor-labelled inhibitory relay; relay inhibits PF_KNEE-F | 0.125 /0.5 |
| TOEDF→PF_FUNC_ANKLE_Z_N | Direct inhibition of the independent dorsiflexor PF channel | 0.68725 /2.749 |

The W2L ancestor put TOEDF inhibition on its shared `PF_KNEE-F`. `functional_gait2392/network.py` explicitly moves that destination to the separate ankle-N channel. That channel's active direct MN output in the **original audit condition only** is `MN_tib_ant` at 0.3 µS; the source mapping selects ankle_df muscles. **The final ground9 output profile instead has TA=0.5161055046 µS and SO=0.03605914573 µS.** The six incoming contact-PF routes and their mappings are unchanged, but this audit is not a measurement or proof of every final-profile output/current. Do not label 0.3 µS as the final ground9 TA gain in a figure or handoff.

The retained toe route **suppresses DF during toe loading**; it does not excite DF or directly suppress knee flexion in this Hunt version. Existing PF_N→its own inhibitory relay remains, so reciprocal effects can also occur. These are graph/function statements, not a predicted net torque under all other inputs.

Heel→InE remains zero in the actual Hunt mask; absent heel→InF routes are not invented. The direct HEEL→contralateral RG-F projection remains6 µS. Restoring the contact-PF slots changes neither of those choices.

The actual caller uses a measured heel event of amplitude0…**2 nA** after preceding whole-foot unloading. Its passive HEEL encoder has a40 ms time constant and a2 mV upper equilibrium. The regional toe port is0.5 times clipped normalized toe load, with an admissible upward-load upper bound **0.75 nA**; the passive TOE encoder likewise has40 ms time constant and0.75 mV upper equilibrium. The actual TOE→TOEDF connection is5 µS at E=8 mV, window0–5 mV; its isolated toe-only upper equilibrium is **3.428571 mV**.

| Static endpoint, recipient at 2.5 mV | Scale0.25 | Scale1 |
| --- | --- | --- |
| HEEL current into each PF_IN, nA | +0.275 | +1.100 |
| Subsequent heel-only relay current into PF_F, nA | −1.570857 | −5.498000 |
| TOEDF current into ankle-DF PF, nA | −3.534429 | −14.137714 |

These endpoints use the actual windows/reversals and isolated relay equilibria. They are not measured delivered currents, global bounds under other PF inputs, or evidence that full-scale contact recruitment is benign. The regional formula has no explicit lower clip, so its nonnegative range assumes upward ground-only contact. Event amplitude is bounded; transient encoder voltage may stay below the equilibrium upper endpoint. Other inputs and relay history must be logged during an actual stage5 test.

The toe pathway can produce strong DF inhibition even at scale 0.25. This warrants checking actual PF_N/tibialis-anterior currents and unloaded clearance during refinement, rather than treating the gain as a harmless graph activation. No gain, event policy, support level, or acceptance gate was changed by this audit.
