# Neuromuscular Control Hypothesis for the Skateboard Kickflip

Working document, 2026-09-08. Companion to `board-sports-gymnastics-motor-control-lit-review.md`.
Target testbeds: **AnimatLab 2** (full neuromechanical sim, BPA muscle plugin), **SNS toolbox** (network architecture/timing standalone), **Simscape Multibody** (plant + optimizer feasibility sweeps).

**Verification tags used throughout** (per Ben, 2026-09-08: be wary of second-hand/partial verification):
- **[V]** — verified peer-reviewed source (see companion lit review / source list below)
- **[L]** — trick-tip lore: multi-source online coaching consensus, NOT peer-reviewed
- **[U]** — unverified estimate (mine or coaching folklore numbers; to be measured, not cited)

---

## 1. The trick, decomposed (with timeline)

| Phase | Duration [U] | Mechanics | Key events |
|---|---|---|---|
| 1. Setup / crouch | 0.2–0.5 s | Countermovement: hip/knee flexion, ankle dorsiflexion; arms back; gaze on board | Stretch-shorten pre-load; postural set |
| 2. Pop / jump | 0.1–0.2 s | Triple extension + rear-ankle plantarflexion slams tail; pop GRF ~2.25 BW [V: Frederick 2006, n=7 elite] | Tail-ground contact breaks; board pitch-up begins |
| 3. Drag → flick | ~50–150 ms after pop | Front foot drags up grip (levels board, as in ollie [V: Nakashima & Chida 2021 — early rapid front-foot pull-up; Wood 2020 — front-knee flexion ↔ board height]), then **ankle** flick off corner of nose [L] imparts roll-axis angular impulse | Board gains ω about its long axis; rider fully unweighted |
| 4. Flight | 0.4–0.7 s | Ballistic COM; board flips ~360° beneath; no rider-board control possible | Predictive phase only; sound/vision cue use plausible [V: Cesari 2014 — sound shapes anticipation/EMG timing] |
| 5. Catch | 0.05–0.15 s | **Back foot catches tail from above** as board levels [L]; front foot held up/out [L]; contact friction + normal force arrest residual roll | Board attitude re-fixed relative to rider |
| 6. Land / absorb | 0.1–0.3 s | Both feet on deck, knee/hip flexion absorbs landing; ollie landing GRF 4.52 ± 0.58 BW at 40–50 ms, forefoot [V: Frederick 2006] | Zero COM velocity without hop (gymnastics "zeroing" framing [V: Pain 2007]) |

Trick-tip consensus (skateboarding.com/Dominick Walker, skatedeluxe, Empire, community threads; all **[L]**): ollie stance with front foot angled just behind front bolts; pop while jumping; **flick with the ankle, not a whole-leg kick**; stay above the board; back foot catches, front foot returns as board levels; the two most-cited failure modes are (a) flicking too hard/low ("kicking the board away") and (b) dropping the back foot early (under-committing the jump).

Back-of-envelope sanity check **[U — my derivation, not citable]**: a standard deck (deck ~0.8 kg, 81 cm × 20 cm, trucks+wheels ~0.7 kg) has roll-axis inertia I_roll ≈ 0.008–0.01 kg·m². A clean kickflip completes ~360° in ~0.4 s → ω ≈ 15 rad/s → required angular impulse L = Iω ≈ 0.15 N·m·s. Over a ~75 ms flick at ~0.4 m moment arm, the tangential flick force is only **order 5 N** — a light, precise, distal impulse, not a powerful kick. This is consistent with the coaching cue "flick with the ankle" and motivates the synergy hypothesis below. (One direct study exists — "Successful KickFlip: Dynamic Analysis..." with an MPU6050 IMU on a miniature board, rigid-body flip model — but it is in a low-tier venue (bcpublication.org); treat as lead, not evidence, until read and sanity-checked.)

---

## 2. Central hypothesis (H0)

> **The kickflip is executed by a discrete, timed feedforward motor program — a sequence of three muscle synergies (pop, drag-flick, catch-land) whose inter-synergy timing is gated by proprioceptive/force afferents at the two brief contact events (tail pop, board catch) — with sensory feedback excluded from the flick itself because pop-to-flick and flip-to-catch intervals are shorter than visual/vestibular loop latency, and with leg impedance set feedforward at landing from an internal model of flight time.**

Rationale from literature:
- Pop→flick→flip→catch each unfold in ~50–150 ms windows; visuomotor loop latency (~150–250 ms) cannot close the loop inside them; and high-gain delayed feedback actively destabilizes board balance (speed wobble from reflex delay [V: Várszegi 2016]; "critical delay" framework [V: Molnár & Insperger 2021/2022]). So the flick must be ballistic/feedforward, with feedback confined to event gating and impedance.
- Expertise signatures elsewhere are predictive feedforward tuning, not more feedback: pre-landing EMG onset 80–90 ms before contact with economical (not maximal) co-activation [V: Pavlasová 2025; Niespodziński 2021; Janshen 2000], and stronger feedforward sensorimotor modulation in trained board athletes [V: Ouyang & Chen 2023].
- The ollie literature already shows strategy/timing, not raw power, drives outcome [V: Heinen 2024 — back-foot-dominant timing +12% height; Nakashima & Chida 2021 — front-foot pull-up timing; Wood 2020 — knee kinematics ↔ height; Candotti 2012 — lower-limb power explains only ~76% of ollie performance].

---

## 3. Sub-hypotheses (each = one ablation/sweep experiment)

### H1 — Afferent-gated timing: the flick fires off the tail-pop afferent, not a free-running clock
**Claim:** Flick onset is triggered/gated by tail-contact force/proprioceptive afferents (short-latency, spinal-ish), making flick timing self-calibrating to pop quality; jitter in pop does not propagate to flip-phase error.
**Test (SNS toolbox → AnimatLab):** Run the network two ways — (a) afferent-gated flick, (b) fixed inter-synergy delay (clock). Sweep pop-magnitude jitter (±20%) and board property perturbations (mass, tail stiffness, friction). 
**Predicts:** (a) maintains flip angle at catch within ±30° over the whole sweep; (b) drifts out of the catch window linearly with pop jitter. 
**Falsified if:** clock-timing matches afferent-gating under all perturbations (then gating adds nothing and H1 is dead).

### H2 — Single-gain flick synergy: the flick is one distal muscle module scaled by one scalar
**Claim:** The flick is a fixed inter-joint synergy (hip flexion + knee extension + ankle plantarflexion/inversion) whose *shape* is constant and whose *amplitude* is set by one gain. Flip rate ω scales ~linearly with flick force × contact duration (angular impulse); contact-point placement sets spin-axis purity (off-corner contact → off-axis/yaw contamination → "varial-ish" or rocket flip outcomes).
**Test (Simscape):** Collapse the front-leg actuation to a synergy vector with scalar gain; sweep gain × contact point on a 2D grid with a proper contact/friction deck model. 
**Predicts:** near-linear ω(gain) map until contact-break saturation; an ellipse-shaped "clean-flip region" in the (gain, contact-point) plane; the known failure taxonomy maps onto the map's boundary (low gain = under-flip, high gain = over-flip/double, medial contact = off-axis spin). 
**Falsified if:** ω is insensitive to gain (contact breaks too early — impulse mechanism wrong) or no clean-flip region exists for ANY fixed synergy shape (then the CNS must shape the synergy per-trial, contradicting the module claim).

### H3 — Active back-foot catch: the catch is a programmed reaching/impedance synergy, not passive waiting
**Claim:** The back leg is programmed (feedforward, flight-time-locked) to extend toward the tail and catch **from above** with tuned impedance; catch contact friction arrests residual roll before landing. The front foot is deliberately held aloft (inhibitory component of the program).
**Test (AnimatLab):** Ablate the back-foot catch synergy (foot stays passive at flight height); and separately ablate the front-foot hold (front foot drops early — the classic beginner error [L]). Measure residual ω at ground contact and landing success over stochastic trials. 
**Predicts:** no-catch ablation collapses success specifically under flick-gain noise (residual ω at landing grows with gain error); early-front-foot ablation produces the "board flips into shins / lands nose-up" failure signature seen in beginners. 
**Falsified if:** passive-both-feet performance is statistically indistinguishable from programmed catch (then the catch needs no neural program).

### H4 — Flight-time-locked stiffness scheduling at landing (internal model + zeroing)
**Claim:** Leg pre-activation for landing is launched from the pop/flight internal model (i.e., timed from pop, corrected only after catch contact), with pre-activation level tuned — economical, not maximal [V: gymnastics landing literature]. 
**Test (AnimatLab/Simscape):** Compare three controllers: (a) pre-activation timed from pop event (flight-time estimate), (b) reactive-only (activation begins at catch/ground contact), (c) (a) + short-latency reflex after contact. Sweep drop height / flight time mismatch (±15%). 
**Predicts:** (a)/(c) keep peak landing GRF near the ~4–5 BW envelope [V] and keep COM-zeroing (no hop); (b) overshoots GRF or requires hop; (a) degrades gracefully with flight-time mismatch until reflex term rescues it. 
**Falsified if:** reactive-only matches predictive timing on all metrics.

### H5 — (Extension) Balance-phase delay ceiling on a rolling setup
**Claim:** For a rolling (non-static) kickflip setup, the pre-pop balance controller has a critical sensory-loop delay below which balance fails, in the structure predicted by the delayed-feedback board-balance literature [V: Várszegi 2016; Molnár & Insperger 2021/2022] — i.e., the feedforward program exists partly to *avoid* operating feedback at high gain near the pop. 
**Test (SNS toolbox, 1-DOF roll model first):** inject loop delay; map stability vs. delay vs. forward speed. **Predicts:** a critical-delay boundary qualitatively matching the Insperger balance-board results. Out of scope for the static flat-ground sim; queue after H1–H4.

### H6 — (Hardware-relevant extension) BPA bandwidth vs. the flick
**Claim:** The 50–100 ms flick timing window is at/below the direct force bandwidth of BPAs; a hardware implementation will need either pressure pre-charge + fast exhaust valving or series elastic release at the ankle — the model should tell us which. 
**Test (Simscape):** swap ideal torque actuators for BPA force models (MonoPam characterization data) in the H2 synergy sweep; check whether any gain achieves the impulse before the board leaves contact. Result feeds actuator/valve requirements for the robot leg.

---

## 4. Testbed plan (platform division of labor)

**Stage 1 — SNS toolbox (network only):** build the program: three synergy modules, their intrinsic timing, the two afferent gates (tail-pop force → flick release; catch contact → landing program), inhibitory front-foot hold. Validate timing logic against the phase table; H1 partially testable here with a 1-DOF board-roll plant.

**Stage 2 — Simscape Multibody (plant + feasibility):** simplified rider (two 3-segment legs + lumped torso, point feet) + deck (rigid body with trucks' mass properties) + unilateral contact with friction (grip tape). Use the Opt_run machinery (surrogateopt → patternsearch) to find the feasibility envelope of synergy parameters — this is "system identification of the motor program": what pop/flick/catch parameters even produce a catchable kickflip. H2 and H6 live here; H4's GRF metrics here.

**Stage 3 — AnimatLab 2 (coupled neuromechanics):** full loop — SNS driving Hill-type or BPA-plugin muscles in the multibody rider+board+ground model; run ablations H1–H4 as the formal hypothesis tests; stochastic trial batches (timing jitter, magnitude noise, board property perturbations) → success-rate statistics.

**Primary metrics:** flip angle at catch (target 360° ± 30°); board attitude (roll/pitch/yaw) at touchdown; residual roll rate at catch (H3); peak landing GRF vs. 4–5 BW envelope (H4); catch success rate over ≥100 stochastic trials per condition; sensitivity ∂(flip error)/∂(flick-delay jitter) (H1).

---

## 5. Parameter anchors (tagged)

| Quantity | Value | Tag | Source |
|---|---|---|---|
| Ollie pop force | ~2.25 BW | V | Frederick et al. 2006, J Appl Biomech |
| Ollie landing peak vGRF | 4.52 ± 0.58 BW @ 40–50 ms, forefoot | V | Frederick et al. 2006 |
| Front-knee flexion ↔ board height | positive correlation | V | Wood et al. 2020, Apunts |
| Back-foot-dominant pop timing | +12% ollie height | V | Heinen et al. 2024, Sports Eng (optimal control, not human experiment) |
| Lower-limb power → ollie performance | ~76% variance | V | Candotti et al. 2012 |
| Pre-landing EMG onset | 80–90 ms before contact (VL) | V | Pavlasová et al. 2025 (scoping review) |
| Reflex-delay-induced board instability | speed wobble = delay in human loop | V | Várszegi et al. 2016, J R Soc Interface |
| Flick force, ω, I_roll, flight/flip times | ~5 N; ~15 rad/s; ~0.01 kg·m²; 0.4–0.7 s / ~0.4 s | U | my back-of-envelope (§1) — must be measured in Stage 2 and/or from video |
| Flick with ankle off nose corner; back-foot catch from above; front foot up | coaching consensus | L | trick tips (below) |
| Back-foot catch first | coaching consensus | L | not peer-reviewed anywhere — H3 is genuinely novel to test |

## 6. Sources for this document

Peer-reviewed (details in companion lit review): Frederick 2006; Wood 2020; Nakashima & Chida 2021; Heinen 2024; Candotti 2012; Cesari 2014; Várszegi 2016; Molnár & Insperger 2021/2022; Pavlasová 2025; Niespodziński 2021; Janshen 2000; Ouyang & Chen 2023; Pain 2007 (abstract-verified).

New this session (use with the stated caution):
- [Low-tier venue — read before trusting] "Successful KickFlip: Dynamic Analysis and Influencing Factors of the Skateboard KickFlip Motion" — IMU-instrumented miniature board, rigid-body flip model. https://bcpublication.org/index.php/SJIRR/article/view/9283 / https://www.researchgate.net/publication/403397798
- [Educational, not peer-reviewed] Tom Rocks Maths / T. Shakespeare, "Exploring Skateboarding Kinematics, World Records and the Intermediate Axis Theorem" (kickflip via angular-momentum conservation; intermediate-axis effects in multi-axis flips). https://tomrocksmaths.com/wp-content/uploads/2024/08/rad-math_-exploring-skateboarding-kinematics-world-records-and-the-intermediate-axis-theorem-thomas-shakespeare.pdf
- [Trick tips, L-tier] Skateboarding.com — How to Kickflip (Dominick Walker): https://www.skateboarding.com/how-to/kickflip-like-a-pro-step-by-step-instruction-tutorial ; Skatedeluxe trick tip: https://www.skatedeluxe.com/blog/en/trick-tips/skateboard/flat/how-to-kickflip/ ; Empire Skate: https://thinkempire.com/blogs/news/empire-trick-tips-comment-faire-un-kick-flip
