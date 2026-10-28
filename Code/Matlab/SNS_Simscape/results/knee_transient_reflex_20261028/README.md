# Transient knee-reflex demonstration

This is an isolated reduced-order SNS demonstration, not a validated human tendon-reflex model. The original configuration does **not** settle: it maintains neural, activation, and mechanical oscillations. The explicitly modified graded-response configuration settles, responds with extension, and returns toward its operating point. No manuscripts, identification coefficients, or biomimetic torque results were changed.

## Run

In MATLAB, add this folder to the path and call `run_transient('replay')`. Call `run_transient('primary')` to repeat the unchanged-parameter diagnostic. Any other tag selects the graded-response comparison. The driver derives paths from its own location, loads the original, and saves only an isolated `KneeTransient.slx` here. It never executes the original destructive builder. `audit_builder` separately runs a snapshot of that builder in `builder_snapshot` and compares against the original loaded model.

`KneeTransient.slx` is runnable after loading its PreLoadFcn parameters; it defaults to the 2 nA pulse. The driver runs the no-stimulus case first and checks settling before running the final pulse experiment. Raw outputs retain the entire 0–50 s startup, baseline, stimulation, and recovery history. Plotting shows 29–40 s with time relative to the 30 s stimulus, without filtering or resetting any state.

## Convention

The original internal equation is positive flexion:

`I*theta_ddot = r*Fflex - r*Fext + Tload - b*theta_dot - K*(theta-th0)`.

At Ben's request, the final figures use the paper convention: `theta_paper=-theta_internal` and `omega_paper=-omega_internal`. Extension is therefore a positive plotted excursion. Raw `log_th`, `log_thd`, and summary arrays C/P/F retain internal coordinates. The JSON includes explicitly named paper-coordinate metrics. This output coordinate transformation does not alter the dynamics.

## Audit and implemented equations

`builder_comparison.txt/.mat`: exact equality of functional dialog parameters and input wiring for every loaded block versus a newly built snapshot, including library internals. `loaded_audit.txt/.mat` captures loaded block parameters and connections. The original library and all six requested inputs are preserved in `original_snapshot`, with SHA256 hashes in `original_sha256.csv`.

Neuron: `dV/dt=1000/Cm*(Gm*(Vrest-V)+Iapp+sum(g*Esyn)-V*sum(g))`, with seconds, mV, nA, nF, and microSiemens. Neural output is `clip((V-Thr)/Slope,0,1)`. Synapses output `[g;g*Esyn]`, with `g=gmax*clip((Vpre-ThrPre)/SlopePre,0,1)`. Despite builder mask labels describing inverse-mV slopes, the actual code divides by these parameters; they act as voltage ramp widths.

Activation: `dA/dt=1000/50*(S-A)`, with output saturation. BPA: `F=clip(Fmax*A*max(0,0.25-strain),0,Fmax)`. Extensor strain is `0.15*theta/ROM`; flexor strain is `0.15*(1-theta/ROM)`. These are simple actuator equations, not the identification model or a physiological muscle model. The NumPy twin was inspected, not executed: it writes to the original result path and uses sequential Euler updates rather than simultaneous continuous integration. Its source and existing result are preserved; it is not used as evidence for the new transient.

## Original-settings failure and explicit changes

The original tonic drives (4.2/2.5 nA), external flexion load (0.5 N m), inertia (0.06 kg m²), damping (0.5 N m s/rad), spring (0.5 N m/rad), rest angle (0.26 rad), ROM (1.5708 rad), moment arm (0.035 m), force scales (500/450 N), and 50 ms activation constant are unchanged throughout.

In the original 25–30 s window, angle range is 0.01242 rad (0.712°), velocity range is 0.6973 rad/s, MN voltage ranges are about 18 mV, activation ranges are about 0.5, and both normalized MN outputs span the entire [0,1] interval. This is persistent switching through saturation, not a settled operating point or simple startup drift. The narrow high-gain saturation ramps are a likely contributor; this experiment does not uniquely identify an instability mechanism. Original pulse results are diagnostic only and cannot establish recovery from rest.

The graded comparison changes exactly six mask parameters: both MN `Slope` values from 1 to 10 mV, and all four Ia synaptic `SlopePre` values from 0.5 to 10 mV. It reduces loop gain while retaining the actual neuron/synapse/activation/force equations, wiring, thresholds, conductances, mechanics, and drives. This is an explicit engineering sensitivity test, not parameter identification or physiological calibration. It was one comparison, not a search for a preferred trace. The final baseline flexor drive is at its zero floor; thus the experiment does not establish reciprocal suppression of a tonically active flexor.

## Protocol and verification

Two time-based steps generate a +2 nA rectangular current from 30.000 to 30.030 s. The pulse is added to the extensor Ia spindle current **after the spindle transducer and before SN_Ia_ext input port 1 (Iapp)**. It then passes through the existing sensory membrane and both Ia synapses to the MNs. It is a tendon-tap-like sensory-current perturbation, not a simulated tendon impact or imposed stretch. Amplitude is a model current, not a measured tendon-tap magnitude.

Both stimulus and no-stimulus cases start from the same original initial conditions and evolve continuously for 50 s. Settling gates use 25–30 s: angle range <1e-5 rad, absolute velocity <1e-5 rad/s, and all six membrane and both activation state ranges <1e-5 in their respective units. Settling is tested on the control before the final stimulus run. Full traces allow independent inspection of earlier settling. A pathway-disabled control was not added; the matched no-stimulus control establishes the incremental stimulus effect without changing baseline feedback.

ODE45 uses relative tolerance 1e-8, absolute tolerance 1e-10, and maximum step 1 ms; the fine run halves the maximum step to 0.5 ms. These are adaptive maximum steps, not fixed timesteps. Raw results contain solver time samples. Quantitative comparisons use linear interpolation to a 1 ms common grid, with no smoothing. Sub-millisecond response timing is not claimed.

The graded baseline is -42.6392° in paper coordinates. Extension peaks at +1.0908° relative to control, 111 ms after pulse onset. The late residual is approximately +9.3e-7° (45–50 s mean versus baseline); the no-stimulus angle range after 25 s is approximately 1.47e-8°. Halving the maximum timestep changes the full angle trace by at most 5.21e-5°. Recovery requires late displacement <0.01° and absolute velocity <1e-4 rad/s; recovery time is the first grid time after the last violation of either bound. See replay_summary.json for the fresh-process result and timing.

Files `primary_*` retain the original-settings failure; `graded_*` retain the first successful comparison; `replay_*` contain the fresh-process verification. Raw MATs include stimulus, MN voltages and normalized outputs, all sensory voltages, activations, forces, angle, velocity, and extensor Ia synaptic pairs. `replay_parameters.mat` captures the final block parameters; console logs contain baseline state ranges. The final editable figure, vector PDF, and PNG are `replay_response.*`.

## Proposed caption

Transient response of the reduced-order SNS knee model to a 2 nA, 30 ms current pulse at the extensor Ia sensory-neuron input after settling. The stimulus case is compared with a matched no-stimulus control. Flexion is negative; the pulse evokes approximately 1.09° of extension, followed by recovery. The displayed graded-response configuration uses broadened Ia synaptic and MN output saturation ramps; the original settings maintained oscillations. This demonstration tests the implemented circuit and actuator equations and does not establish physiological reflex fidelity.
