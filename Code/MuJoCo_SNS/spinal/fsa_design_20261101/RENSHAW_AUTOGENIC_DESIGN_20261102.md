# Separate autogenic Renshaw candidate and local transfer proof

Session 1, easteregg2. This note records the initial isolated design of **Hunt_PT_RC_autogenic_0p1mV_v1**, a separately identified optional neural-wiring candidate. Its later default-off integration and independent fresh stages 1–4 plus quarter-strength stage-5 body acceptance are indexed in [RENSHAW_INTEGRATION_REVIEW_20261102.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/RENSHAW_INTEGRATION_REVIEW_20261102.md); those results were not part of the isolated proof below. Historical constructors/results remain preserved. Tendons, damping, geometry, and other mechanics are unchanged.

## Observed omission and explicit choice

The selected trial23 `physical_eval.build` reconstructs 1,033 actual cells with the exact ordered names in its saved mechanics. Four extra entries in the name dictionary are logging aliases. All92 inherited RC→MN pairs exist but have zero live conductance after motor masking; none is restored by a saved PF/bi/MN/reflex key. In contrast, 92 RC→IaIN and 4,140 RC→RC outputs remain active at 0.5 µS. The generic motor-input clearing policy, not the absence of Renshaw cells, suppresses direct recurrent inhibition. No separate RC-specific biological rationale was found. [Selected-caller audit](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/input_policy_v1/HUNT_RENSHAW_DISPOSITION.json>).

The new candidate restores only the existing **autogenic motor-collateral feedback** slots, with an explicit modest isolated engineering contract. It does not blindly restore the inherited 0.5 µS or identify a physiological Renshaw gain. Motor collateral is not an afferent; its independent buffer remains active in the explicitly deafferented stage-1 protocol.

## Derivation using actual retained parameters

Let \(u\) be MN voltage and \(v\) RC voltage, measured in mV relative to rest. The actual source pair has \(G_M=G_R=1\) µS, zero rest/bias, collateral \(h=1\) µS with \(E_C=8\) mV, feedback \(E_R=-5\) mV, and both activation windows 0–5 mV. Define \(a(u)=\operatorname{clip}(u/5,0,1)\) and \(b(v)=\operatorname{clip}(v/5,0,1)\). The isolated recurrent operator is

\[
C_M\dot u=-G_Mu+I+h_R b(v)(E_R-u),\qquad
C_R\dot v=-G_Rv+h a(u)(E_C-v).
\]

Here \(h_R\) is the new feedback conductance; membrane capacitance is in µF, current in nA, and time in seconds. The actual source capacitances are 0.030 µF for RC and either 0.030 or 0.045 µF for MN. They are preserved, not accelerated.

Apply baseline \(I=G_Mu_0=2.5\) nA, giving \(u_0=2.5\) mV without feedback. The declared **closed-loop** target is \(u_*=2.4\) mV, a −0.10 mV change. Therefore

\[
v_*={h a(u_*)E_C\over G_R+h a(u_*)}
={0.48\cdot8\over1.48}=2.5945945946\ \mathrm{mV},
\]
\[
h_R={G_M(u_0-u_*)\over b(v_*)(u_*-E_R)}
=0.0260416666667\ \mu\mathrm{S}.
\]

At that equilibrium the feedback current is exactly −0.10 nA. This derivation includes conductance shunting and the collateral's dependence on the resulting MN voltage; it does not hold RC artificially at its no-feedback voltage. All92 anatomically paired operators derive the same feedback conductance. Other incoming circuits are silent **only in this unit-contract definition** and remain intact in the full candidate.

## Local state-space and transfer function

Within the interior activation branch, at \((u_*,v_*)\), the Jacobian is

\[
A=\begin{bmatrix}
-(G_M+h_R b_*)/C_M & h_R b'_*(E_R-u_*)/C_M\\
h a'_*(E_C-v_*)/C_R & -(G_R+h a_*)/C_R
\end{bmatrix},\quad
B=\begin{bmatrix}1/C_M\\0\end{bmatrix},\quad
C=\begin{bmatrix}1&0\end{bmatrix},\ D=0.
\]

The driving-force derivatives are retained. The small-signal transfer from applied MN current to MN voltage is

\[
H(s)=C(sI-A)^{-1}B
={B_1(s-A_{22})\over(s-A_{11})(s-A_{22})-A_{12}A_{21}}.
\]

For the 0.030 µF MN case,

\[
A\approx\begin{bmatrix}-33.783784&-1.284722\\36.036036&-49.333333\end{bmatrix}\ {s^{-1}},
\quad H(s)\approx{33.333333s+1644.444444\over s^2+83.117117s+1712.962963}.
\]

| Actual MN capacitance | Continuous poles, s⁻¹ | Actual 2 ms Euler poles | DC gain, mV/nA |
| --- | --- | --- | --- |
| 0.030 µF | −37.796800, −45.320317 | 0.9244064, 0.9093594 | 0.960000 |
| 0.045 µF | −23.727899, −48.127956 | 0.9525442, 0.9037441 | 0.960000 |

These are the **isolated two-cell operator's** poles and transfer, not full-network poles, 105-coordinate body stability, contact stability, or a gait Floquet analysis. Simultaneous RC lateral inhibition and other full-network inputs can change the operating point and effective feedback substantially.

## Completed electrical proof

The actual-source two-cell SNS fixture has freely evolving MN and RC voltages and a constant applied current; no voltage clamp or new population is introduced into the full controller. Both actual capacitance cases converge within 5.8×10⁻¹⁵ mV of the declared equilibrium. Native dense versus accelerated sparse output is exact. The analytic sampled Jacobian \(I+0.002A\) matches finite differences of the actual SNS update within 2.43×10⁻¹¹. A 1 pA, 100 ms input pulse agrees with the linearized prediction within 7.28×10⁻⁸ mV.

On the actual selected full neural constructor:

- Only the92 existing autogenic RC output conductances change. All other neuronal fields, retained RC→IaIN/RC→RC routes, and existing learned gain buffers are exact.
- Commands 0/0.5/1/0 preserve those outputs through a **separate motor-collateral buffer**, independently of afferent masks. Other command-dependent routes match the unchanged comparator exactly.
- Existing integrated local/timing/coincidence checks remain valid. Stage-1 deafferentation still verifies while RC collateral persists.
- The full 1,033-cell network matches native dense versus accelerated sparse stepping exactly for 128 electrical updates.
- Deliberate sparse, dense-cache, and authoritative-buffer corruptions are rejected.
- An independent fresh process reproduces both serialized JSON files byte-for-byte.

The full-network stepping test is electrical only. It does not establish whole-network phase effects, physiological gain identification, improved knee behavior, or body acceptance.

## Artifacts and opt-in use

Implementation: [renshaw_autogenic_design.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/renshaw_autogenic_design.py). Proof runner: [verify_renshaw_autogenic.py](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/verify_renshaw_autogenic.py).

Exact candidate [PROFILE.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/renshaw_autogenic_v1/PROFILE.json>), [PROOF.json](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/renshaw_autogenic_v1/PROOF.json>), and [fresh proof](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/integration/renshaw_autogenic_v1/fresh/PROOF.json>). Profile SHA256: `10e944bba8455ebd9277b1aa0c2857216f9504a677d95c4336d55fff38fdb867`; proof SHA256: `19b114efe9315f426b08a4ab34e29ec1c1d75987abad5c82624f274366f9b43b`.

After constructing the selected existing integrated neural candidate, call `renshaw_autogenic_design.install(net, profile_path)` before body binding. `install(net, None)` is an exact no-op. `verify(net)` checks its separate authoritative conductances/cache. The installer rejects an existing unknown instance-level command wrapper, unpaired source parameters, non-passive MN/RC mechanisms, an incompatible source operator, or already active feedback slots. It adds no input, neuron, force, trajectory, or clamp. Following the initial isolated proof and explicit parental release, a default-off caller hook was added; no controller/balance/physics parameter changed. [Exact integration and coactive full 184-state qualification](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/RENSHAW_INTEGRATION_REVIEW_20261102.md). In the complete coactive MN/RC subsystem, retained45-peer inhibition reduces the isolated−0.10 mV effect to approximately−0.00630 mV. Its full 184-state sampled poles are locally stable, without implying whole-controller/body stability.

Before any adoption, bind/calibrate the real PT object with this new neural mapping, then recheck the same air/standing/ground criteria and declared support forces as a **distinct candidate**. No previous body pass is transferred, and no new body or training job was launched here.
